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

otherwise try to install the packages with the following commands:
```
sudo apt install [package_name] (like dssp, mafft, fasttree)
```
On MacOS:
```
brew install [package_name] (like brewsci/bio/dssp, mafft, fasttree)
```

# Nucleotide Structure Database (nustruDB)
## Complete data construction pipeline
To construct and filter the data in a more fast forward way, a pipeline script is provided. But alternatively if only parts of the pipeline are needed, the scripts can be run individually as described in the following sections. Several improvement have been made to parallelize and asynchronize the data fetching over API calls. Even pagination was used to limit the request as applicable. However, many website like ncbi limit API requests to around 3-10 requests per second and a account is needed. Dealing with hundred thousands of uniprot entries, this can slow down the process.

Note: Some of the script require an api key for ncbi to accelerate and parallelize the data fetching. The api key can be created following the website instructions: https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us. The scripts are optimized to use the api key and the mail address to fetch the data, if requested changes can be made to lower the request rate and avoid the api key.

The database construction follows the following steps for individual entries:
- Fetch entry from PDB or Uniprot
- Map the nucleotide sequence to the protein sequence (Uniprot is seperated in two steps)
- Filter the entries for wrong translations and redundancy
- Assign secondary structure and other features to the entries

The database construction follows the following steps for protein families (Only for Uniprot):
- Fetch protein members of family from InterPro
- Fetch entries from Uniprot
- rest is same as for individual entries

While the individual scripts are described in the following sections, there is a complete pipeline script called `nustrufiller.py`. Based on the input single or list of uniprot id(s) or intepro id(s), it will fetch and construct the required data. The pipelin script can be run with the following command:
```
usage: nustrufiller.py [-h] -i INPUT -o OUTPUT_PATH [-n NAME] [-u UNIQUE] [-w] -m API_MAIL -k API_KEY [-d]

Complete pipeline to fetch the data for uniprot IDs or interpro IDs (complete protein families).

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Single uniprotID or interproID, or file with comma seperated ids.
  -o OUTPUT_PATH, --output OUTPUT_PATH
                        Output path or directory to store the log file and the data.
  -n NAME, --name NAME  Name of the output files and log file.
  -u UNIQUE, --unique UNIQUE
                        pssibility to drop duplicate. To keep duplicates use None. Default: organisms
  -w, --overwrite       If file name already exists, overwrite it. Default is False.
  -m API_MAIL, --mail API_MAIL
                        Provide mail for ncbi api account: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/.
  -k API_KEY, --key API_KEY
                        Provide api key from ncbi api account: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/.
  -d, --download        Download the structure files.
```

Note: This pipeline script only works for Uniprot data. For PDB data, the script in "Create nustruDB for PDB data" section should be used.

## Individual steps to create nustruDB database 
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

##### Slow version (but complete)
````
usage: UPmapNT [-h] -i ENTRYID [--sql] [--pandas] -o OUTPUT_PATH -n NAME [--create-fasta {protein,nucleotide,all}] [-w] -m API_MAIL -k API_KEY

Map uniprot ID to nucleotide sequence and prints an allignment of the pdb protein sequence to the nucleotide sequence

options:
  -h, --help            show this help message and exit
  -i ENTRYID, --input ENTRYID
                        Path to file containing PDB IDs with comma separation.
  --sql                 Store the data in a SQL database.
  --pandas              Store the data in a pandas DataFrame.
  -o OUTPUT_PATH, --output OUTPUT_PATH
                        Output path or directory to store the log file and the data.
  -n NAME, --name NAME  Name of the output files and log file.
  --create-fasta {protein,nucleotide,all}
                        Create a fasta file with the nucleotide sequences, protein sequences or both.
  -w, --overwrite       If file name already exists, overwrite it. Default is False.
  -m API_MAIL, --mail API_MAIL
                        Provide mail from ncbi api account: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
  -k API_KEY, --key API_KEY
                        Provide api key from ncbi api account: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
````
Test this uniprot mapping with the following command:
````
 python fastUPmapNT.py -i Example/examples_nustruDB/example_uniprotList.tsv -o . -n example_uniprotList_nustru [-w]
````
##### Fast version (but needs additional steps) - recommended
To use the parallel version of the uniprot mapping, a tsv table of the uniprot ids and the column feature is needed. The script is adapted from: https://www.uniprot.org/help/api_queries and uses batches to retrieve the required data. The script is called with the following command:
```
usage: fetchUPTSV.py [-h] -i INPUT_FILE -o OUTPUT_PATH -n NAME [-w]

Retrieve uniprot features from uniprot IDs and store them in a tsv file

options:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input INPUT_FILE
                        Input file with uniprot IDs.
  -o OUTPUT_PATH, --output OUTPUT_PATH
                        Output path or directory to store the data.
  -n NAME, --name NAME  Name of the output file.
  -w, --overwrite       If file name already exists, overwrite it. Default is False.
```
Test the tsv table creation with the following command:
```
 python fetchUPTSV.py -i Example/examples_nustruDB/example_uniprotIDs.txt -o uniprot_features -n uniprot_features [-w]
```

The tsv table can then be used to map the nucleotide sequences to the uniprot ids with the following command:
````
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
````

Test the uniprot mapping with the following command:
```
python fastUPmapNT.py -i Example/examples_nustruDB/example_uniprotList.tsv -o . -n example_uniprotList_nustru [-w]
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
<br />

### Assign secondary structure and other features to the entries (both PDB and Uniprot)
It is better to filter first the data and then fetch the secondary structure information, as it reduced the API load. To run this script, an installed version of the `DSSP` software is required. Try to install it with the following command.
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
python db_fetch.py -i Example/examples_nustruDB/example_dbfetch.csv -o . -n example_dbfetch_nustru_fetched [-w]
```
<br />

## Other scripts to obtain data
### Fetch the protein members of a protein family from InterPro
Some of the analysis was done on whole protein families. While `nustrufiller.py` can be used to fetch the data and perform several data preperation steps automatically, the script `fetchInterPro.py` can be used to fetch the protein members of a protein family from InterPro family seperatelöy. The script is called with the following command:

```
usage: fetchINPRO.py [-h] -o OUTPUT_PATH [-w] family_id

Fetch Accessions from InterPro

positional arguments:
  family_id             InterPro family ID to fetch accessions for

options:
  -h, --help            show this help message and exit
  -o OUTPUT_PATH, --output OUTPUT_PATH
                        Output to store the interpro accessions.
  -w, --overwrite       If file name already exists, overwrite it. Default is False.
```
Test the interpro fetching with the following command:
```
python fetchINPRO.py IPR000839 -o IPR000839_accessions.txt
```
<br />
<br />

# Codon and protein analysis
## Correlation and cross-validation analysis of codon usage bias and protein structure
To analyse the effect of codon usage bias and specific correlation between synonymous codons and secondary structure elements ...

## Protein family based analysis of codon rarity, secondary structure and evolution
The evolutionary analysis is based on the multiple sequence alignment of different protein families with data from the nustruDB (described previously). For each MSA the amino acid normalized codon rarities are calculated with the following formula:
$$CRS_{column}= {\sum \limits _{AA_{AA}} ^{n_{aln}}(\sum \limits _{occ=1} ^{n_{aln}} {AA_{occ}} * f_c) \over {len_{total}(alignment)} - (gaps)}$$ 
<br />

Where $AA_{AA}$ is the amino acid at the position of the alignment, $AA_{occ}$ is the amino acid at the position of the alignment, $f_c$ is the frequency of the codon that encodes the amino acid, $len_{total}(alignment)$ is the total length of the alignment and $gaps$ are the gaps in the alignment. The script is called with the following command:

### Make the multiple sequence alignment

The multiple sequence alignment is done with the `mafft` software. The software can be installed with the following command:
On Linux:
```
sudo apt install [package_name] (like dssp, mafft, fasttree)
```
On MacOS:
```
brew install [package_name] (like brewsci/bio/dssp, mafft, fasttree)
```

The previous database scripts should have created a protein fasta file with the filtered protein sequences, so that the alignment can then be done like this:
```
mafft --auto input.fasta > output.fasta
```

### Analysis of codon rarities, structure and evolutionary implications
For the analysis a jupyter notebook `/protein_family_analysis/msa_codon_evol.ipynb`is provided. The notebook explains each individual step and the results of the analysis. It directly implements the described formula to calculate the codon rarity on the msa. The notebook is divided into the following sections:
- Calculate and map the Codon rarity at each position of a multiple sequence alignment
- Transform alignment and sequences to arrays and align the gaps if necessary
- Calculate the codon rarity for each position of the alignment
- Visualize the Codon Rarity Score for each position of the alignment in MSA
- Calculate the correlation between the smoothed Y and the secondary structure frequencies as a log odds ratio (rare codons and secondary structure)
