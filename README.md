# BiologicalData_project

# PF03060 Protein Family Analysis Pipeline

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

## Folder Tree

<pre> Biological_Data_project/
│
├── 01_model_building.ipynb
├── 02_ground_truth_and_characterization.ipynb
├── notebook_01_pipeline_old.ipynb
│
├── data/
│ ├── raw_PF03060_uniprot.tsv
│ ├── PF03060_full_length.fasta
│ ├── PF03060_metadata.csv
│ ├── swissprot_pos_PF03060.tsv
│ ├── swissprot_neg_PF03060.tsv
│ └── ...
│
├── results/
│ ├── Q12723_10-372_clustalo_msa.fasta
│ ├── Q12723_10-372_clustalo_msa_NR90.fasta
│ ├── Q12723_10-372_clustalo_msa_NR90_clean.fasta
│ ├── PF03060_PSSM_logodds.csv
│ ├── PF03060_profileHMM.json
│ ├── PF03060_MSA_stats.csv
│ ├── PF03060_eval_scores.csv
│ └── ...
│
└── README.md
``` </pre>

# Notebook 1 – Model Construction  
`01_model_building.ipynb`

This notebook implements the first major block of the project.

## Steps Performed

1. Retrieve homologous sequences
   - BLASTP search of seed domain (Q12723, positions 10–372)
   - UniRef50 database via EBI REST API

2. Build Multiple Sequence Alignment
   - Clustal Omega
   - 201 sequences aligned

3. Redundancy Reduction
   - 90% identity threshold (NR90)
   - Reduced to 152 sequences

4. Clean MSA
   - Remove columns with >50% gap fraction
   - Alignment reduced from 1620 columns to 350 columns

5. Build Statistical Models
   - PSSM (log-odds vs background frequencies)
   - Lightweight profile HMM (match states + emission probabilities)

6. Conservation Analysis
   - Gap profile
   - Shannon entropy per column
   - Identification of highly conserved residues

Outcome:
- Clean domain alignment (~350 residues)
- Reproducible statistical models saved in `results/`

---

# Notebook 2 – Model Evaluation and Characterization  
`02_ground_truth_and_characterization.ipynb`

This notebook evaluates the predictive performance of the models.

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
Step 1: Fetch UniProt sequences for PF03060
        Extract Pfam domain sequences
        Build MSA, clean it, generate PSSM and HMM
notebooks/notebook_01_pipeline



### Output
<pre> ```
data/ — input/output sequences and metadata
results/ — MSAs, cleaned MSAs, PSSM, HMM files
``` </pre>
