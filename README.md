# BiologicalData_project

# PF03060 Protein Family Analysis Pipeline

This project builds a workflow to analyze the Pfam family **PF03060**, from fetching sequences to building MSAs, PSSMs, and HMMs.

## Folder Tree

<pre> ``` Biological_Data_project/ │ ├── notebook_01_pipeline.ipynb ├── README.md │ ├── data/ │ ├── raw_PF03060_uniprot.tsv │ ├── PF03060_full_length.fasta │ ├── PF03060_metadata.csv │ ├── pf03060_domains.fasta │ └── domain_metadata.csv │ ├── results/ │ ├── pf03060_msa_raw.fasta │ ├── pf03060_msa_clean.fasta │ └── ... │ └── scripts/ └── utils.py ``` </pre>


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
data/ — input/output sequences and metadata
results/ — MSAs, cleaned MSAs, PSSM, HMM files

