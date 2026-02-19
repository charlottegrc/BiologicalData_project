from __future__ import annotations
import csv
import time
from dataclasses import dataclass
from collections import Counter, defaultdict, deque
from typing import Dict, Iterable, List, Optional, Tuple, Set
import requests
import xml.etree.ElementTree as ET
from src.session_utils import EBI_TAXONOMY, NS, UNIPROT_REST, make_session
import networkx as nx
import matplotlib.pyplot as plt


def load_accessions(accessions: Iterable[str]) -> List[str]:
    """
    Accepts a Python list (or any iterable) of UniProt accessions.
    
    Example:
        load_accessions(['O43099', 'P56577', 'P56578'])
    """
    if accessions is None:
        raise ValueError("You must provide a list of accessions.")

    if isinstance(accessions, str):
        raise TypeError("Provide a list of accessions, not a single string.")

    accs = [str(a).strip() for a in accessions if a and str(a).strip()]
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


def prune_to_labeled_and_ancestors(G: nx.DiGraph, node_counts, label_top_n: int) -> nx.DiGraph:
    """Keep only top-N nodes (by count) and all their ancestors to preserve hierarchy."""
    top = {n for n, _ in node_counts.most_common(label_top_n)}
    keep = set()

    for n in top:
        if n in G:
            keep.add(n)
            keep.update(nx.ancestors(G, n))  # parents/grandparents... up to root

    return G.subgraph(keep).copy()

def plot_taxonomy_hierarchical_pruned(
    node_counts,
    edges,
    out_path="taxonomy_tree_hier_pruned.png",
    max_nodes=200,      # first filter: only consider top max_nodes in the universe
    label_top_n=50,     # second filter: show only these as labeled "targets"
    fig_size=(20, 12),
    font_size=8
):
    # 1) restrict universe to top max_nodes by count
    keep = {n for n, _ in node_counts.most_common(max_nodes)}
    filtered_edges = [(u, v) for (u, v) in edges if u in keep and v in keep]

    # 2) build graph from edges only (prevents isolates)
    G = nx.DiGraph()
    G.add_edges_from(filtered_edges)

    if G.number_of_nodes() == 0:
        raise ValueError("No nodes to plot. Increase max_nodes or check your edges.")

    # 3) find roots, create super-root if needed
    has_parent = {v for (_, v) in G.edges()}
    roots = sorted([n for n in G.nodes() if n not in has_parent])

    if len(roots) == 1:
        root = roots[0]
    else:
        root = "__ROOT__"
        G.add_node(root)
        for r in roots:
            G.add_edge(root, r)

    # 4) prune graph: keep only top label nodes + their ancestors
    G = prune_to_labeled_and_ancestors(G, node_counts, label_top_n=label_top_n)

    # If pruning removed everything (rare), stop cleanly
    if G.number_of_nodes() == 0:
        raise ValueError("Pruning removed all nodes. Try increasing label_top_n/max_nodes.")

    # 5) recompute reachable from root (in case pruning disconnected things)
    reachable = set()
    q = deque([root]) if root in G else deque([n for n in G.nodes() if G.in_degree(n) == 0])
    while q:
        u = q.popleft()
        if u in reachable:
            continue
        reachable.add(u)
        for v in G.successors(u):
            q.append(v)
    G = G.subgraph(reachable).copy()

    # 6) compute depth for hierarchical coordinates (BFS)
    # pick a root that exists
    if root not in G:
        root = sorted([n for n in G.nodes() if G.in_degree(n) == 0])[0]

    depth = {root: 0}
    q = deque([root])
    while q:
        u = q.popleft()
        for v in G.successors(u):
            if v not in depth:
                depth[v] = depth[u] + 1
                q.append(v)

    levels = defaultdict(list)
    for n, d in depth.items():
        levels[d].append(n)
    max_depth = max(levels.keys())

    # 7) positions: y = -depth, x spread in each layer
    pos = {}
    for d in range(max_depth + 1):
        level_nodes = sorted(levels[d], key=lambda n: node_counts.get(n, 0), reverse=True)
        k = len(level_nodes)
        if k == 1:
            xs = [0.0]
        else:
            mid = (k - 1) / 2.0
            xs = [i - mid for i in range(k)]

        spacing = 1.8
        for x, n in zip(xs, level_nodes):
            pos[n] = (x * spacing, -d)

    # 8) sizes + labels (now everything left is meaningful)
    sizes = []
    for n in G.nodes():
        c = node_counts.get(n, 0)
        sizes.append(max(120, min(2500, c * 250)))

    plt.figure(figsize=fig_size)
    nx.draw_networkx_edges(G, pos, arrows=False, alpha=0.25, width=1.0)
    nx.draw_networkx_nodes(G, pos, node_size=sizes, alpha=0.85)

    # Label all nodes in pruned graph EXCEPT the artificial root
    labels = {n: f"{n}\n({node_counts.get(n, 0)})" for n in G.nodes() if n != "__ROOT__"}

    nx.draw_networkx_labels(
        G, pos, labels=labels, font_size=font_size,
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7)
    )

    plt.axis("off")
    plt.tight_layout()
    plt.savefig(out_path, dpi=250)
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

    write_lineage_tsv(lineages, f"{out_prefix}_lineage.tsv")
    write_node_counts_tsv(node_counts, f"{out_prefix}_node_counts.tsv")

    # plot (optional)
    # if node_counts and edges and nx is not None and plt is not None:
    #     plot_taxonomy_graph(node_counts, edges, out_path=f"results/taxonomy/{out_prefix}_tree.png")
    return node_counts, edges