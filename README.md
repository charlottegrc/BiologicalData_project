# BiologicalData_project

# PF03060 Protein Family Analysis Pipeline

This project builds a workflow to analyze the Pfam family **PF03060**, from fetching sequences to building MSAs, PSSMs, and HMMs.

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

### Usage
Step 1: Fetch UniProt sequences for PF03060
        Extract Pfam domain sequences
        Build MSA, clean it, generate PSSM and HMM
notebooks/notebook_01_pipeline



### Output
data/ — input/output sequences and metadata
results/ — MSAs, cleaned MSAs, PSSM, HMM files

