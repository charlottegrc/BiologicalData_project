from __future__ import annotations
import csv
import time
from dataclasses import dataclass
from collections import Counter
from typing import Dict, Iterable, List, Optional, Tuple, Set
import requests
import xml.etree.ElementTree as ET

from src.session_utils import EBI_TAXONOMY, NS, UNIPROT_REST, make_session
try:
    import networkx as nx
    import matplotlib.pyplot as plt
except Exception:
    nx = None
    plt = None

def load_accessions(
    accessions: Optional[Iterable[str]] = None,
    input_path: Optional[str] = None,
    column_candidates: Tuple[str, ...] = ("accession", "uniprot_acc", "uniprot_id", "id"),
) -> List[str]:
    """
    Accepts either:
      - a list/iterable of accessions
      - a TSV/CSV file containing an accession column (any of column_candidates)
      - a plain text file with one accession per line
    """
    if accessions is not None:
        accs = [a.strip() for a in accessions if a and str(a).strip()]
        return sorted(set(accs))

    if not input_path:
        raise ValueError("Provide either 'accessions' or 'input_path'")

    # Try TSV/CSV first
    if input_path.lower().endswith((".tsv", ".csv")):
        delim = "\t" if input_path.lower().endswith(".tsv") else ","
        with open(input_path, "r", newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter=delim)
            if not reader.fieldnames:
                raise ValueError("Empty TSV/CSV file")

            col = None
            lower_map = {c.lower(): c for c in reader.fieldnames}
            for cand in column_candidates:
                if cand.lower() in lower_map:
                    col = lower_map[cand.lower()]
                    break
            if col is None:
                raise ValueError(f"Could not find accession column in {reader.fieldnames}")

            accs = []
            for row in reader:
                v = (row.get(col) or "").strip()
                if v:
                    accs.append(v)
            return sorted(set(accs))

    # Otherwise treat as plain list file
    with open(input_path, "r", encoding="utf-8") as f:
        accs = [line.strip() for line in f if line.strip() and not line.startswith("#")]
    return sorted(set(accs))


# ---------- Lineage extraction ----------

@dataclass(frozen=True)
class TaxonNode:
    name: str

# Downloads the UniProt XML file for a protein.
# ex: https://rest.uniprot.org/uniprot/P12345.xml
# raw XML text (string)
# The project asks: extract lineage from entry/organism/lineage in UniProt XML
def fetch_uniprot_xml(session: requests.Session, accession: str) -> str:
    # same endpoint style as api_and_parsing.py: /uniprot/{ACC}.xml
    url = f"{UNIPROT_REST}/uniprot/{accession}.xml"
    r = session.get(url)
    r.raise_for_status()
    return r.text


# Extracts the lineage terms from XML.
# This gives us the ordered lineage hierarchy
def parse_lineage_from_uniprot_xml(xml_text: str) -> List[TaxonNode]:
    """
    Project requirement: entity/organism/lineage in UniProt XML.
    In UniProt XML, lineage terms typically appear as <lineage><taxon>...</taxon>...</lineage>
    We return ordered lineage nodes (root -> ... -> leaf-ish).
    """
    root = ET.fromstring(xml_text)

    # Try the canonical lineage path
    taxa = root.findall(
        "uniprot:entry/uniprot:organism/uniprot:lineage/uniprot:taxon",
        NS
    )
    if taxa:
        return [TaxonNode(t.text.strip()) for t in taxa if t is not None and t.text]

    # Fallback: sometimes XML variants differ; try any lineage/taxon anywhere under organism
    taxa = root.findall(".//uniprot:organism//uniprot:lineage//uniprot:taxon", NS)
    return [TaxonNode(t.text.strip()) for t in taxa if t is not None and t.text]

# Extracts the NCBI taxonomy ID.
# If XML lineage is missing or incomplete, we can use this ID to fetch lineage from EBI API
def parse_taxon_id_from_uniprot_xml(xml_text: str) -> Optional[str]:
    """
    Optional helper: get NCBI taxonomy id if present.
    UniProt XML often has:
      <dbReference type="NCBI Taxonomy" id="9606"/>
    """
    root = ET.fromstring(xml_text)
    for ele in root.findall("uniprot:entry/uniprot:organism/uniprot:dbReference", NS):
        if ele.attrib.get("type") == "NCBI Taxonomy" and ele.attrib.get("id"):
            return ele.attrib["id"]
    return None


# Backup method if XML parsing fails
def fetch_lineage_from_ebi(session: requests.Session, taxon_id: str) -> List[TaxonNode]:
    """
    Uses EBI Proteins API taxonomy lineage endpoint (as in api_and_parsing.py).
    Returns ordered nodes as provided by the endpoint.
    """
    r = session.get(f"{EBI_TAXONOMY}/lineage/{taxon_id}", headers={"Accept": "application/json"})
    r.raise_for_status()
    data = r.json()
    nodes = []
    for node in data.get("taxonomies", []):
        nm = node.get("scientificName") or node.get("commonName")
        if nm:
            nodes.append(TaxonNode(str(nm)))
    return nodes


# Orchestrates the full lineage retrieval.
def get_lineage(
    session: requests.Session,
    accession: str,
    prefer_xml: bool = True,
    allow_ebi_fallback: bool = True,
    polite_sleep: float = 0.0,
) -> List[TaxonNode]:
    """
    Main function to obtain lineage for one protein accession.
    By default, it uses UniProt XML (project requirement) and can fallback to EBI lineage if needed.
    """
    if polite_sleep:
        time.sleep(polite_sleep)

    xml_text = fetch_uniprot_xml(session, accession)

    if prefer_xml:
        lineage = parse_lineage_from_uniprot_xml(xml_text)
        if lineage:
            return lineage

    if allow_ebi_fallback:
        taxid = parse_taxon_id_from_uniprot_xml(xml_text)
        if taxid:
            lineage2 = fetch_lineage_from_ebi(session, taxid)
            if lineage2:
                return lineage2

    return []


# ---------- Build abundance tree ----------
# Count how many proteins include each taxon
# Create parent-child relationships
def build_tree_edges(lineages: Dict[str, List[TaxonNode]]) -> Tuple[Counter, Set[Tuple[str, str]]]:
    """
    Returns:
      - node_counts: how many proteins include this node in their lineage
      - edges: parent->child edges aggregated across all proteins
    """
    node_counts = Counter()
    edges: Set[Tuple[str, str]] = set()

    for acc, lin in lineages.items():
        # count nodes once per protein (avoid overweighting repeated nodes)
        seen = set()
        for node in lin:
            if node.name not in seen:
                node_counts[node.name] += 1
                seen.add(node.name)

        # edges along the lineage path
        for parent, child in zip(lin[:-1], lin[1:]):
            edges.add((parent.name, child.name))

    return node_counts, edges


# Saves detailed lineage per protein.
def write_lineage_tsv(lineages: Dict[str, List[TaxonNode]], out_path: str) -> None:
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["accession", "lineage_index", "taxon_name"])
        for acc, lin in lineages.items():
            for i, node in enumerate(lin):
                w.writerow([acc, i, node.name])

# Saves abundance per taxon.
def write_node_counts_tsv(node_counts: Counter, out_path: str) -> None:
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["taxon_name", "count"])
        for name, c in node_counts.most_common():
            w.writerow([name, c])



def plot_taxonomy_graph(
    node_counts,
    edges,
    out_path="taxonomy_tree.png",
    max_nodes=200,
):
    

    # Keep only most frequent nodes (avoid huge unreadable graph)
    keep = {n for n, _ in node_counts.most_common(max_nodes)}
    filtered_edges = [(u, v) for (u, v) in edges if u in keep and v in keep]

    G = nx.DiGraph()
    G.add_nodes_from(list(keep))
    G.add_edges_from(filtered_edges)

    # Node size proportional to abundance
    sizes = [max(200, node_counts[n] * 200) for n in G.nodes()]

    pos = nx.spring_layout(G, seed=42)

    plt.figure(figsize=(14, 10))

    # Draw edges
    nx.draw_networkx_edges(G, pos, arrows=False, alpha=0.3)

    # Draw nodes
    nx.draw_networkx_nodes(G, pos, node_size=sizes, alpha=0.8)

    # Create labels including counts
    labels = {n: f"{n}\n({node_counts[n]})" for n in G.nodes()}

    nx.draw_networkx_labels(G, pos, labels=labels, font_size=8)

    plt.axis("off")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.show()



# ---------- Main pipeline ----------

def run_taxonomy(
    accessions: List[str],
    out_prefix: str = "taxonomy",
    polite_sleep: float = 0.0,
) -> None:
    session = make_session()

    lineages: Dict[str, List[TaxonNode]] = {}
    for acc in accessions:
        try:
            lin = get_lineage(session, acc, prefer_xml=True, allow_ebi_fallback=True, polite_sleep=polite_sleep)
            if lin:
                lineages[acc] = lin
        except Exception as e:
            # keep going; log minimal info
            print(f"[WARN] {acc}: {e}")

    node_counts, edges = build_tree_edges(lineages)

    write_lineage_tsv(lineages, f"data/taxonomy/{out_prefix}_lineage.tsv")
    write_node_counts_tsv(node_counts, f"data/taxonomy/{out_prefix}_node_counts.tsv")

    # plot (optional)
    if node_counts and edges and nx is not None and plt is not None:
        plot_taxonomy_graph(node_counts, edges, out_path=f"results/taxonomy/{out_prefix}_tree.png")
    return node_counts, edges