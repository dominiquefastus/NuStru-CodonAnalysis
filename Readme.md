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
Both PDB and Uniprot provide coding sequences for the deployed proteins, but follow different mapping strategies to the nucleotide sequence.

To create the database from PDB, a graphql api is used to fetch the nucleotide sequences. The coding sequence of the protein is annotated with the start of the nucleotide sequence, the end of the nucleotide sequence, the strand and eventual exon shifts from the NCBI nt database. The database is created with the following command:

```
usage: PDBmapNT [-h] -i ENTRYID [--sql] [--pandas] -o OUTPUT_PATH -n NAME [--map-uniprot]

Map PDB ID to nucleotide sequence and prints an allignment of the pdb protein sequence to the nucleotide sequence

options:
  -h, --help            show this help message and exit
  -i ENTRYID, --input ENTRYID
                        Path to file containing PDB IDs with comma separation.
  --sql                 Store the data in a SQL database.
  --pandas              Store the data in a pandas DataFrame.
  -o OUTPUT_PATH, --output OUTPUT_PATH
                        Output path or directory to store the log file and the data.
  -n NAME, --name NAME  Name of the output files and log file.
  --map-uniprot         Map uniprot ID to nucleotide sequence.
```

Test the pdb mapping with the following command:
```
python PDBmapNT.py -i nustruDB/Example/example1_pdbIDs.txt --pandas -o /home/usr/ -n pdb_example1
```
\
To create the database from Uniprot, the nucleotide sequence is fetched from the NCBI nt database. The refered Genebank ID or EMBL ID is used to fetch the nucleotide sequence. Two scripts are provided, but while the first one follows a similiar strategy as the pdb mapping with only providing uniprot ids (see `nustruDB/Example/example1_uniprotIDs.txt`), the second script takes in a predefined list (see `nustruDB/Example/example1_uniprotList.tsv`). The list has to be downloaded manualy at the uniprot database page. The database is then created with the following command:

```
usage: 2fastUPmapNT.py [-h] -i INPUT_FILE -o OUTPUT_PATH -n NAME

Retrieve nucleotide sequences from uniprot IDs.

options:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input INPUT_FILE
                        Input file with uniprot IDs.
  -o OUTPUT_PATH, --output OUTPUT_PATH
                        Output path or directory to store the log file and the data.
  -n NAME, --name NAME  Name of the output files and log file.
  ```

Test the uniprot mapping with the following command:
````
python 2fastUPmapNT.py -i nustruDB/Example/example1_uniprotList.tsv -o /home/usr/ -n uniprot_example1
````

## Evolutionary Analysis of Protein and Nucleotide Sequences
