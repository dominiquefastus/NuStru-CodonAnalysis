# NuStru - Nucleotide Codon Usage Bias Analysis on Protein Structure & Evolution

Mainly, the project aims to discover correlations between specific motifs in the protein structure and a bias in the codons that encode these motifs. It also explores the evolutionary implications of these biases on protein families. 

## Requirements
The project is developed in Python 3.11 or later. All packages are installed respictively in a conda environment. An environment with the required packages is available in the `nustru_environment.yml` file. To create the environment, the following command can be run:
```
conda env create -f setup/nustru_environment.yml
```
For the phylogenetic tree construction and analysis `ete3` is used, which has some specific installation dependencies. So the package need to be installed in an extra environment. The environment with the required packages is available in the `nustru-phyl_environment.yml` file. To create the environment, the following command can be run:
```
conda env create -f setup/nustru-phyl_environment.yml
```
One additional package needs to be installed manually, as the CAI python package has some problems with the pip installation. The package can be installed with the following command:
```
pip install git+https://github.com/Benjamin-Lee/CodonAdaptationIndex.git
```
<br />
Some of the scripts require also external software packages to be installed. Which includes:

```
dssp                    4.4.0
FastTree                2.1.11
mad                     2.2 
mafft                   7.525
```

While the installation procedure for these packages is mentioned in the specific sections, there is also a script provided to install all the required software packages at once. The script can be run with the following command:
```
bash setup/package_helper.sh 
```

otherwise try to install the packages with the following commands:
On Linux:
```
sudo apt install [package_name] (like dssp, mafft, fasttree, mad)
```
On MacOS:
```
brew install [package_name] (like brewsci/bio/dssp, mafft, fasttree, mad)
```
<br />

Sometimes mad can not be installed with these methods, so the source code can be downloaded from the website and installed manually. The source code can be found here: https://www.mikrobio.uni-kiel.de/de/ag-dagan/ressourcen
<br />
<br />

# Nucleotide Structure Database (nustruDB)
## Complete data construction pipeline
To construct and filter the data in a more fast forward way, a pipeline script is provided. But alternatively if only parts of the pipeline are needed, the scripts can be run individually as described in the following sections. Several improvements have been made to parallelize and asynchronize the data fetching over API calls. Even pagination was used to limit the requests if applicable. However, many websites like NCBI limit API requests to around 3-10 requests per second. Dealing with hundreds of thousands of uniprot entries, this can slow down the process.

Note: Some of the scripts require an api key for NCBI to accelerate and parallelize the data fetching. The api key can be created by following the website instructions: https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us. The scripts are optimized to use the api key and the mail address to fetch the data, if requested changes can be made to lower the request rate and avoid the api key.

The database construction includes the following steps for individual entries:
- Fetch entry from PDB or Uniprot
- Map the nucleotide sequence to the protein sequence (Uniprot is seperated in two steps)
- Filter the entries for wrong translations and redundancy
- Assign secondary structure and other features to the entries

The database construction includes the following steps for protein families (only for Uniprot):
- Fetch protein members of a family from InterPro
- Fetch entries from Uniprot
- other steps are the same as for individual entries

While the individual scripts are described in the following sections, there is a complete pipeline script called `nustrufiller.py`. Based on the input, one or a list of Uniprot id(s) or Intepro id(s), it will fetch and construct the required data. The pipeline script can be run with the following arguments / options:
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

Note: This pipeline script only works for Uniprot data. For PDB data, the script in the "Create nustruDB for PDB data" section should be used.

## Individual steps to create nustruDB database 
Both PDB and Uniprot provide coding sequences for the deployed protein structures, but follow different mapping strategies to the nucleotide sequence.
<br />

#### Create nustruDB for PDB data
To create the database from PDB, a graphql api is used to fetch the nucleotide sequences. The coding sequence of the protein is annotated with the start of the nucleotide sequence, the end of the nucleotide sequence, the strand and eventual exon shifts from the NCBI nt database. The csv file is created with the following arguments / options:

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
To create the database from Uniprot, the nucleotide sequence is fetched from the NCBI nt database. The refered Genebank ID or EMBL ID is used to fetch the nucleotide sequence. Two scripts are provided, but while the first one follows a similiar strategy as the PDB mapping with only providing Uniprot ids (see `nustruDB/Example/example1_uniprotIDs.txt`), the second script takes in a predefined list (see `nustruDB/Example/example1_uniprotList.tsv`).
The database is then created with the following arguments / options:

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
To use the parallel version of the uniprot mapping, a tsv table of the uniprot ids and the column features is needed. The script is adapted from: https://www.uniprot.org/help/api_queries and uses batches to retrieve the required data. 
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

The tsv table can then be used to map the nucleotide sequences to the uniprot ids with the following arguments / options:
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
The dataframe or database contains many redundant entries. To reduce the redundancy and also delete wrong translations from the nucleotide coding sequence, the entries are filtered by the following criteria: The coding sequence starts with ATG, the coding sequence is dividable by trinucleotides (codons), the coding sequence only contains valid nucleotides (A, T, C, G) and the coding sequence can be translated to the same assigned protein sequence. The script is called with the following arguments / options:

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
It is better to filter first the data and then fetch the secondary structure information, as it reduces the API load. To run this script, an installed version of the `DSSP` software is required. Try to install it with the following command:
On Linux:
```
sudo apt install dssp
```
On MacOS:
```
brew install brewsci/bio/dssp
```
<br />

The script will create two additional columns with the b-factor for PDB entries or the pLDDT score for Uniprot entries  and the secondary structure. The chain of the model is taken into account. The script is called with the following arguments / options:

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

### Fetch the domains and fold classes
The domains and fold classes were not fetched during the previous steps. The script `db_fclass.py` is used to fetch the annotated domains and fold classes for a protein. The fold class is retrieved from
the CATH or SUPFAM database, this still needs to be improved. The script is called with the following arguments / options:
```
usage: db_fclass.py [-h] -i INPUT_FILE -o OUTPUT_PATH -n NAME [-w]

Fetch the annotated domains and fold classes for a protein.

options:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input INPUT_FILE
                        Input file of csv formatted uniprot or pdb entries.
  -o OUTPUT_PATH, --output OUTPUT_PATH
                        Output to store the new csv with domain and fold class information.
  -n NAME, --name NAME  Name of the output files and log file.
  -w, --overwrite       If file name already exists, overwrite it. Default is False.
```
Test the domain and fold class fetching with the following command:
```
python -i Example/examples_nustruDB/example_dbfclass.csv -o . -n example_dbfclass_nustru_fetched
```
<br />

## Other scripts to obtain data
### Fetch the protein members of a protein family from InterPro
Some of the analysis were done on whole protein families. While `nustrufiller.py` can be used to fetch the data and perform several data preperation steps automatically, the script `fetchINPRO.py` can be used to fetch the members (Uniprot IDs) of a protein family from InterPro seperately. The script is called with the following arguments / options:

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
Test the InterPro fetching with the following command:
```
python fetchINPRO.py IPR000839 -o IPR000839_accessions.txt
```
<br />
<br />

# Codon and protein analysis
## Correlation and cross-validation analysis of codon usage bias and protein structure
To analyse the effects of codon usage bias and specific correlations between synonymous codons and secondary structure elements, two notebooks are provided. The first notebook `codon_metrics_analysis.ipynb` investigates the secondary structure and codon usage bias of different parts of the protein sequence. 
The second notebook `cross_validation_destribution_analysis.ipynb` uses the KL divergence to compare the synonymous codon usage bias of the secondary structure elements. The notebook also compares the frequencies and probabilities. 

## Protein family based analysis of codon rarity, secondary structure and evolution
The evolutionary analysis is based on the multiple sequence alignment of different protein families with data from the nustruDB (described previously). For each msa, the amino acid`s normalized codon rarities are calculated with the following formula:
<br />
<p align="center">
  <img src="https://latex.codecogs.com/svg.image?\inline&space;\large&space;\bg{white}$$CRS_{column}={\sum\limits&space;_{AA_{AA}}^{n_{aln}}(\sum\limits&space;_{occ=1}^{n_{aln}}{AA_{occ}}*f_c)\over{len_{total}(alignment)}-(gaps)}$$" title="$$CRS_{column}={\sum\limits _{AA_{AA}}^{n_{aln}}(\sum\limits _{occ=1}^{n_{aln}}{AA_{occ}}*f_c)\over{len_{total}(alignment)}-(gaps)}$$" />
</p>
<br />

Where $AA_{AA}$ is the amino acid at the position of the alignment, $AA_{occ}$ is the amino acid at the position of the alignment, $f_c$ is the frequency of the codon that encodes the amino acid, $len_{total}(alignment)$ is the total length of the alignment and $gaps$ are the gaps in the alignment.

### Make the multiple sequence alignment

The multiple sequence alignment is done with the `mafft` software. The software can be installed with the following command:
On Linux:
```
sudo apt mafft
```
On MacOS:
```
brew install mafft
```

The previous database scripts should have created a protein fasta file with the filtered protein sequences, so that the alignment can then be done like this:
```
mafft --auto input.fasta > output_aln.fasta
```

### Tree construction and rooting
The tree construction is done with the `FastTree` software. The software can be installed with the following command:
On Linux:
```
sudo apt install FastTree
```
On MacOS:
```
brew install FastTree
```
From the previous alignment, the tree can be constructed with the like this:
```
FastTree input_aln.fasta > output_tree.nwk
```
The automatic rooting of the tree can be done with MAD. The software should be manually installed from the source code, as described in the requirements section. The rooting can be done with the following command:
```
python mad.py input_tree.nwk 
```

Alternatively, the msa and tree construction can be automatically done with the `setup/align_and_tree.sh` script or the `setup/alignment_and_manual_tree_rooting.ipynb` notebook. Both helps to align the sequences and construct the tree. The notebook also helps to root the tree manually if desired.

### Analysis of codon rarities, structure and evolutionary implications
For the analysis of the codon rarity and protein structure the notebook `msa_codon_structure.ipynb` is provided. The notebook explains each individual step of the analysis. It directly implements the described formula to calculate the codon rarity on the msa. Furthermore, the log odds ratio is calculated to reveal positions with a higher likelihood to be influenced by codon usage bias. 
<br />
<br />
For the phylogenetic analysis the notebook `phyl_codon_rarity.ipynb`is provided. The notebook explains the individual steps and the results of the analysis. It uses the `ete3` package to visualize the tree and the codon rarity. It also does a regression analysis to find the correlation between the codon rarity and the evolutionary distance.

## Protein fold analysis of domains and fold classes with codon usage bias
To analyse the fold classes and domains on the codon rarity (CRS) or relative synonymous codon usage (RSCU) the notebook `fold_domain_analysis.ipynb` is used. The notebook visualizes the ratio of alpha helix and beta sheet contents and the CRS or RSCU. To investigate the CRS on the folds, the path to the pickle file from `msa_codon_structure.ipynb`, which should be created when running the notebook, needs to be defined.

# Implications and bugs
Some of the scripts and notebooks are planned to be more optimized. The logging for the data construction seems to not work properly and should be fixed.