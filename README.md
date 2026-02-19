# BiologicalData_project

This project is designed to run in a Jupyter Notebook environment, like Google Colab or amother Linux-based Jupyter setup. The notebook installs and executes external bioinformatics tools (e.g., Clustal Omega) using shell commands (apt-get), which may not work in standard Python scripts or non-Linux systems

# PF08534 Protein Family Analysis Pipeline

# Project Overview

The goal of this project is to:

1. Build statistical models of the PF03060 domain
   - Position-Specific Scoring Matrix (PSSM)
   - Lightweight profile Hidden Markov Model (HMM)

2. Evaluate these models against SwissProt proteins
   - Compare predictions to official Pfam annotations
   - Compute classification metrics (Precision, Recall, F1, MCC)

3. Characterize the biological properties of the domain family
   - Taxonomic distribution
   - Functional annotation (GO terms)
   - Conservation analysis


# Notebook 01_02 – Model Construction and Evaluation and Characterization  
`01_02_Construction_Evaluation.ipynb`


## Ground Truth Definition

SwissProt proteins were divided into:

- Positives:
  - Reviewed proteins annotated with Pfam PF03060
- Negatives:
  - Reviewed proteins without PF03060 annotation

This provides a controlled evaluation dataset.


## Requirements

### Python Packages
- Python 3.9+
- Biopython
- Pandas
- NumPy

Install with:

```bash
pip install biopython pandas numpy
```

### External Tools

Clustal Omega (MSA) — http://www.clustal.org/omega/
HMMER (hmmbuild for HMM) — http://hmmer.org/
Optional: JalView (visualize MSA) — http://www.jalview.org/
Optional: BLAST+ (homology search) — ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

## Usage
Step 1: Fetch UniProt sequences for PF08534
        Extract Pfam domain sequences
        Build MSA, clean it, generate PSSM and HMM
notebooks/notebook_01_pipeline

### Output
<pre> ```
data/ — input/output sequences and metadata
results/ — MSAs, cleaned MSAs, PSSM, HMM files
``` </pre>

**Notebook 01_02 - Construction & Evaluation**

- **File:** [01_02_Construction_Evaluation.ipynb](01_02_Construction_Evaluation.ipynb)
- **Purpose:** Build MSA-based models (PSSM and HMM) from BLAST/UniRef hits, clean alignments, generate models, run searches vs SwissProt, and compute protein- and residue-level evaluation metrics against Pfam ground truth.

- **Key data inputs:**
   - `data/O43099_Uniref50_1000.xml` (BLAST XML of UniRef50 hits)
   - SwissProt FASTA for evaluation (`uniprot_sprot.fasta`)

- **Main steps (high level):**
   - Parse BLAST XML and map pairwise alignments to query-length sequences.
   - Write a combined FASTA with the query first and mapped subject sequences.
   - Build an MSA with Clustal Omega (keeps query first with `--output-order=input-order`).
   - Clean the MSA: normalize characters, remove gappy sequences/columns, deduplicate, and write `O43099_final_msa.fasta`.
   - Build a profile HMM with `hmmbuild` and a PSSM with `psiblast -in_msa ... -out_pssm`.
   - Create an ungapped BLAST DB and run `psiblast` and `hmmsearch` against SwissProt to collect predictions.
   - Parse `hmm` domtblout and `pssm` TSV outputs to produce per-protein and per-residue prediction sets.
   - Compute confusion matrices and metrics (Precision, Recall, F1, Balanced Accuracy, MCC) at both protein and residue (domain) levels.

- **Commands used (examples from the notebook):**

   - Install tools (Colab / Debian):

      ```bash
      apt-get -qq update && apt-get -qq install -y clustalo hmmer ncbi-blast+
      pip install biopython pandas
      ```

   - Run Clustal Omega to produce MSA (keeps query first):

      ```bash
      clustalo -i data/O43099_Uniref50_1000.fasta -o data/O43099_raw_msa.fasta --force --outfmt=fasta --output-order=input-order
      ```

   - Build HMM and search SwissProt:

      ```bash
      hmmbuild data/O43099_Uniref50_1000.hmm data/O43099_final_msa.fasta
      hmmsearch --tblout data/O43099_hmm_vs_swissprot.tbl --domtblout data/O43099_hmm_vs_swissprot.domtbl data/O43099_Uniref50_1000.hmm uniprot_sprot.fasta
      ```

   - Build PSSM and run psiblast against SwissProt:

      ```bash
      # makeblastdb for an ungapped FASTA
      makeblastdb -in data/protFamily_db_ungapped.fasta -dbtype prot -out data/blastdb/protFamily_db_ungapped

      psiblast -in_msa data/O43099_final_msa.fasta -db data/blastdb/protFamily_db_ungapped -num_iterations 1 -out_pssm data/O43099_Uniref50_1000.pssm -out data/O43099_pssm.log
      psiblast -in_pssm data/O43099_Uniref50_1000.pssm -db uniprot_sprot.fasta -outfmt "6 sseqid qstart qend sstart send evalue bitscore length pident" -out data/O43099_pssm_vs_swissprot.tsv
      ```

- **Outputs produced by the notebook:**
   - Cleaned MSA: `data/O43099_final_msa.fasta`
   - Profile HMM: `data/O43099_Uniref50_1000.hmm`
   - PSSM: `data/O43099_Uniref50_1000.pssm`
   - Search results: `data/O43099_hmm_vs_swissprot.domtbl`, `data/O43099_pssm_vs_swissprot.tsv`
   - Evaluation tables and metrics printed in the notebook and exported to `results/` as CSVs (confusion matrices, per-model metrics).

- **Notes & reproducibility tips:**
   - The notebook supports running in Google Colab (mount Drive) or locally; set `data_dir = 'data'` for local runs.
   - The cleaning step enforces an allowed amino-acid alphabet, upper-cases sequences, pads/truncates to equal length, removes sequences/columns by gap fraction, and deduplicates.
   - When running large database searches (SwissProt), ensure sufficient local disk and memory; consider running on a machine with BLAST/HMMER installed or using Colab with small test sets.
