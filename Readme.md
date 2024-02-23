# NuStru - Nucleotide Codon Usage Bias Analysis on Protein Folding & Structure

Mainly, the project aims to discover correlations between specific 3D motifs in the protein structure and a bias in the condons that encode these motiffs. It also explores the evolutionary implications of these biases. 

## Requirements
The project is developed in Python 3.11 or later. The following packages are required to run the project:
```
python                    3.11.6 
biopython                 1.81
mysql-connector-python    8.3.0
numpy                     1.26.0 
pandas                    2.2.0 
pandarallel               1.6.5
tqdm                      4.66.1
```

All packages are installed respictively in a conda environment. An environment with the required packages is provided in the `nustru-environment.yml` file. To create the environment, the following command can be run:
```
conda env create -f nustru-environment.yml
```

## Setup to create nustruDB database
Both PDB and Uniprot provide coding sequences for the deployed proteins, but follow different mapping to the nucleotide sequence.
