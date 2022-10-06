# Subepitope Tools

Subepitope Tools is a Python library for analyzing antibody epitopes and paratopes with a focus on analysis of the various features contained within. 

## Usage

Copy ```config_example.ini``` to ```config.ini``` in the configuration folder. In your ```config.ini``` file set your project path and data path.

This repository was setup to automate the analysis of structural information present in antibody-antigen complexes which are deposited on [RCSB PDB](https://www.https://www.rcsb.org/). 

### General Guidelines
1. Determine how to define your antigen of interest. For example, it could be a uniprot ID or sequence. Please follow the [RCSB PDB Search API](https://search.rcsb.org) guidelines on how to define your antigen definition and search. 

### Usage as a Library

modules can be used for custom analysis. For example, if a user simply wants to pull a report from the PDB containing chain information of several PDBs of interest, they can use functions to classify those chains as heavy chain, light chain or light chain.  

```python
from subepitope-tools import preprocess_structures

# Appends heavy chain, light chain or antigen to RCSB summary csv
df['peptide_label'] = df.apply(preprocess_structures.peptide_label_df, axis = 1)
```
