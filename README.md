# BiologicalData_project

# PF03060 Protein Family Analysis Pipeline

This project builds a workflow to analyze the Pfam family **PF03060**, from fetching sequences to building MSAs, PSSMs, and HMMs.

## Folder Tree

<pre> ``` Biological_Data_project/
│
├── notebook_01_pipeline.ipynb      # Your main pipeline notebook (Step 1–5)
├── README.md                       # Instructions, description, software used
│
├── data/                           # Raw and processed data files
│   ├── raw_PF03060_uniprot.tsv     # Full UniProt entries (Step 1)
│   ├── PF03060_full_length.fasta   # Extracted full-length sequences (Step 2)
│   ├── PF03060_metadata.csv        # Metadata for sequences (Step 2)
│   ├── pf03060_domains.fasta       # Domain sequences (Step 2)
│   └── domain_metadata.csv         # Domain metadata (Step 2)
│
├── results/                        # Results from analyses (MSA, HMM, etc.)
│   ├── pf03060_msa_raw.fasta       # Crude seed-based MSA (Step 3)
│   ├── pf03060_msa_clean.fasta     # Cleaned MSA (Step 3b)
│   └── ...                         # Additional analysis outputs (e.g., HMMs)
│
└── scripts/                        # Optional: separate Python helper scripts
    └── utils.py                     # Functions for fetching, processing sequences
``` </pre>


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
