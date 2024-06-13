# NuStru - Nucleotide Codon Usage Bias Analysis on Protein Folding & Structure

Mainly, the project aims to discover correlations between specific motifs in the protein structure and a bias in the codons that encode these motifs. It also explores the evolutionary implications of these biases. 

## Requirements
The project is developed in Python 3.11 or later. The following packages are required to run the scripts in the project:
```
python                    3.11.6 
biopython                 1.81
biopandas                 0.5.1.dev0 
mysql-connector-python    8.3.0
numpy                     1.26.0 
pandas                    2.2.0 
pandarallel               1.6.5
tqdm                      4.66.1
```

All packages are installed respictively in a conda environment. An environment with the required packages is provided in the `nustru-environment.yml` file. To create the environment, the following command can be run:
```
conda env create -f setup/nustru-environment.yml
```
<br />
Some of the scripts require also external software packages to be installed. Which indludes:

```
dssp                    4.4.0
mmseqs2                 15-6f452
mafft                   7.525
FastTree                2.1.11
```

While the installation procedure for these packages is mentioned in the specific sections, there is also a script provided to install all the required software packages at once. The script is run with the following command:
```
bash setup/package_helper.sh 
```

## Complete data construction pipeline
To construct and filter the data in a more fast forward way, the following pipeline is provided. But alternatively if only parts of the pipeline are needed, the scripts can be run individually as described in the following sections.



## Setup to create nustruDB database
Both PDB and Uniprot provide coding sequences for the deployed proteins, but follow different mapping strategies to the nucleotide sequence.
<br />

#### Create nustruDB for PDB data
To create the database from PDB, a graphql api is used to fetch the nucleotide sequences. The coding sequence of the protein is annotated with the start of the nucleotide sequence, the end of the nucleotide sequence, the strand and eventual exon shifts from the NCBI nt database. The csv file is created with the following command:

```
usage: PDBmapNT [-h] -i ENTRYID [--sql] [--pandas] -o OUTPUT_PATH -n NAME [--map-uniprot] [--create-fasta {protein,nucleotide,all}] [-w]

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
  --create-fasta {protein,nucleotide,all}
                        Create a fasta file with the nucleotide sequences, protein sequences or both.
  -w, --overwrite       If file name already exists, overwrite it. Default is False.
```

Test the pdb mapping with the following command:
```
python PDBmapNT.py -i Example/examples_nustruDB/example_pdbIDs.txt --pandas -o . -n example_pdbIDs_nustru --create-fasta all
``` 

#### Create nustruDB for Uniprot data
To create the database from Uniprot, the nucleotide sequence is fetched from the NCBI nt database. The refered Genebank ID or EMBL ID is used to fetch the nucleotide sequence. Two scripts are provided, but while the first one follows a similiar strategy as the pdb mapping with only providing uniprot ids (see `nustruDB/Example/example1_uniprotIDs.txt`), the second script takes in a predefined list (see `nustruDB/Example/example1_uniprotList.tsv`).
The database is then created with the following command:

```
usage: fastUPmapNT.py [-h] -i INPUT_FILE -o OUTPUT_PATH -n NAME [-w]

Retrieve nucleotide sequences from uniprot IDs in a csv file.

options:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input INPUT_FILE
                        Input file with uniprot IDs in a csv format.
  -o OUTPUT_PATH, --output OUTPUT_PATH
                        Output path or directory to store the log file and the data.
  -n NAME, --name NAME  Name of the output files and log file.
  -w, --overwrite       If file name already exists, overwrite it. Default is False.
```

Test the uniprot mapping with the following command:
```
python fastUPmapNT.py -i Example/examples_nustruDB/example_uniprotList.tsv -o . -n example_uniprotList_nustru [-w]
```
<br />

### Assign secondary structure and other features to the entries (both PDB and Uniprot)
To run this script, an installed version of the `DSSP` software is required. Try to install it with the following command.
On Linux:
```
sudo apt install dssp
```
On MacOS:
```
brew install brewsci/bio/dssp
```
<br />

The script will create two additional columns with the b-factor for pdb entries or the pLDDT score for uniprot entries as dictionary by position and the secondary structure as a 1-dimensional sequence. The chain of the model is taken into account. The script is called with the following command:

```
usage: db_fetch.py [-h] -i INPUT_FILE -o OUTPUT_PATH -n NAME [-d] [-w]

Fetch the pdb and map the secondary structure, as well as plddt or bfactor to the protein sequence.

options:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input INPUT_FILE
                        Input file of csv formatted uniprot or pdb entries.
  -o OUTPUT_PATH, --output OUTPUT_PATH
                        Output to store the new csv with secondary structure information.
  -n NAME, --name NAME  Name of the output files and log file.
  -d, --download        Download the pdb files.
  -w, --overwrite       If file name already exists, overwrite it. Default is False.
```
Test the secondary structure mapping with the following command:
```
python db_fetch.py -i Example/examples_nustruDB/example_dbfetch.csv -o . -n example_dbfetch_nustru_fetched
```
<br />


### Filter the entries for the analysis
The dataframe or database contains many redundant entries. To reduce the redundancy and also wrong translations from the nucleotide coding sequence, the entries are filtered by the following criteria. The coding sequence start with ATG, the coding sequence is dividable by a trinucleotide (codon), the coding sequence only contains valid nucleotides (A, T, C, G) and the coding sequence can be translated to the same assigned protein sequence. The script is called with the following command:

```
usage: db_filter.py [-h] -i INPUT_FILE -o OUTPUT_PATH -n NAME [-u UNIQUE] [-w]

Filter the csv databases to reduce rendundancy, check cds and right translation and further improve the data contained.

options:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input INPUT_FILE
                        Input file of csv formatted uniprot and/ or pdb entries.
  -o OUTPUT_PATH, --output OUTPUT_PATH
                        Output to store the new filtered and reduced csv.
  -n NAME, --name NAME  Name of the output csv file and log file.
  -u UNIQUE, --unique UNIQUE
                        pssibility to drop duplicate. To keep duplicates use None. Default: organisms
  -w, --overwrite       If file name already exists, overwrite it. Default is False.
```
Test the filtering with the following command:
```
python db_filter.py -i Example/examples_nustruDB/example_db_filtered.csv -o . -n nustruDB -u 
```

## Evolutionary Analysis of Protein and Nucleotide Sequences

$$ CR_{position}= {\sum \limits _{AA_{AA}} ^{n_{aln}}(\sum \limits _{occ=1} ^{n_{aln}} {AA_{occ}} * f_c) \over {len_{total}(alignment)} - (gaps)} $$