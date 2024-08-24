#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
INFO: Pandarallel will run on 10 workers.
INFO: Pandarallel will use standard multiprocessing data transfer (pipe) to transfer data between the main process and workers.
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
  
  Script to fetch the annotated domains and fold classes for a protein from the UniProt database. The fold classes are assigned based on the SCOP and CATH databases, however this is still under development. 
  The script requires the fetchUPTSV.py script to fetch the UniProt data, it reads a csv file with the protein sequences and the UniProt IDs. 
  Example python -i Example/examples_nustruDB/example_dbfclass.csv -o . -n example_dbfclass_nustru_fetched
  """

import re
import shutil
import argparse

import logging
import pandas as pd

from pathlib import Path

from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True, nb_workers=10, verbose=0) 

from fetchUPTSV import UniProtFetcher

class sequenceDomainFetcher:
    def __init__(self):
        self.domain_pattern = re.compile(r'DOMAIN (\d+)\.\.(\d+); /note="([^"]+)"; /evidence="([^"]+)"')
        self.scop_fold_classes = {"1000000":"All alpha proteins", "1000001":"All beta proteins", "1000002":"Alpha and beta proteins (a/b)", 
                                  "1000003":"Alpha and beta proteins (a+b)", "1000004":"Small proteins"}
        self.cath_fold_classes = {"1":"Mainly Alpha", "2":"Mainly Beta", "3":"Alpha Beta", "4":"Few Secondary Structures", "6":"Special"}
    
    def extract_domains(self, text):
        """Extract the domain information from the column"""
        # Regular expression to match the domain information
        domain_pattern = self.domain_pattern
        
        # Find all matches in the text
        matches = domain_pattern.findall(text)
        
        # Create a list of lists for each domain
        domains = []
        for match in matches:
            start, end, note, evidence = match
            domains.append([start, end, note, evidence])
        
        return domains

    def assign_foldclass(self, data, output_path, name):
        """Assign the fold class to the protein"""
        # this method is not complete yet and will be extended
        # first try to assign domains that are annotated in the SCOP
        # hit = scopdb[scopdb["SF-UNIID"] == f"{data["primary_id"]}"]["SF-UNIREG"]

        # extract the domains from the domain column
        domains = self.extract_domains(data['domain']) 
        
        # the gene3d database contains direct information about the fold class
        # split the gene3d column by the semicolon and remove the quotes
        # get a list of lists with the gene3d information
        gene3d = data['gene3d'].replace('"','').split(';')
        gene3d = [gene3d[i:i + 3] for i in range(0, len(gene3d), 3)]
        
        # the supfam database contains direct information about the fold class
        # however the database is not retievable with API requests nor is it downloadable
        # supfam = data['supfam'].replace('"','').split(';')
        # supfam = [supfam[i:i + 3] for i in range(0, len(supfam), 3)]
        
        # loop through the domains and get the fold class for each domain
        if len(domains) > 0 and len(domains[0]) > 0: 
            cath_fclass = []
            # supfam_fclass = []      
            
            # only get the fold class for the first domain
            # this needs to be changed to get the fold class for all domains            
            if len(gene3d[0]) > 0:
                for gene in gene3d:
                    if len(gene) > 1:
                        cath_fclass.append(self.cath_fold_classes.get(gene[0][0], "Unknown"))
            
            """if len(supfam[0]) > 0:
                for sup in supfam:
                    if len(sup) > 1:
                        supfam_fclass.append(sup[0])"""
                        
            family = data['family'].replace('.','').replace('"','').split(';')
            family.remove('') if '' in family else family
            family = [name for name in family if not name.startswith('IPR')]


            # write the new data to the csv file and delete the data and new data to free memory
            new_data = pd.DataFrame([{'source': data['source'], 'primary_id': data['primary_id'], 'gene_name': data['gene_name'], 'organism': data['organism'], 'expression_system': data['expression_system'],
                    'protein_sequence': data['protein_sequence'], 'nucleotide_id': data['nucleotide_id'], 'nucleotide_sequence': data['nucleotide_sequence'], 'bfactor_or_plddt': data['bfactor_or_plddt'],
                    'secondary_structure': data['secondary_structure'], 'family': family, 'domains': domains, 'cath_fclass': list(set(cath_fclass))}])
            new_data.to_csv(f'{output_path}/{name}.csv', mode='a', index=False, header=False)
            del data, new_data
        
def main():
    parser = argparse.ArgumentParser(
        prog='db_fclass.py',
        description="Fetch the annotated domains and fold classes for a protein."
    )
    parser.add_argument(
        '-i', '--input', type=str, dest="input_file", required=True,
        help='Input file of csv formatted uniprot or pdb entries.'
    )
    parser.add_argument( 
        '-o', '--output', type=str, dest="output_path", required=True,
        help='Output to store the new csv with domain and fold class information.'
    )
    parser.add_argument(
        '-n', '--name', type=str, dest="name", required=True,
        help='Name of the output files and log file.'
    )
    parser.add_argument(
        '-w', '--overwrite', action="store_true", dest="overwrite", required=False, default=False,
        help='If file name already exists, overwrite it. Default is False.' 
    )
    args = parser.parse_args()
    
    logging.basicConfig(filename=f'{args.output_path}/{args.name}.log',
            filemode='a',
            format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
            datefmt='%H:%M:%S',
            level=logging.ERROR)
    
    # read the csv file as input
    nucleotide_protein_seqs_df = pd.read_csv(args.input_file)

    # create the new csv file with the respective columns if it does not exist or the overwrite flag is set
    if not Path(f'{args.output_path}/{args.name}.csv').exists() or args.overwrite:
        with open(f'{args.output_path}/{args.name}.csv', mode='w') as f:
            f.write('source,primary_id,gene_name,organism,expression_system,protein_sequence,nucleotide_id,nucleotide_sequence,bfactor_or_plddt,secondary_structure,family,domains,cath_fclass\n')
    else:
        print(f"Error: {args.output_path} already exists. Use -w to overwrite.")
        exit(1)

    uniprotIDs = nucleotide_protein_seqs_df["primary_id"].values.tolist()

    UPfetcher = UniProtFetcher()
        # submit the ID mapping job and get the job ID
    job_id = UPfetcher.submit_id_mapping(from_db="UniProtKB_AC-ID", to_db="UniProtKB", ids=uniprotIDs)

    # check if the ID mapping results are ready and get the results
    if UPfetcher.check_id_mapping_results_ready(job_id):
        link = UPfetcher.get_id_mapping_results_link(job_id)
        # get the ID mapping results from the search URL and store them in a tsv file with the defined features as columns
        # %2Cxref_supfam_full
        results = UPfetcher.get_id_mapping_results_search(link + "?format=tsv" + "&fields=accession%2Cxref_pdb_full%2Cft_domain%2Cxref_gene3d_full%2Cxref_interpro_full") 
        
    # store the results in a tsv file with the defined features as columns 
    columns = results[0].split("\t") 
    data = [line.split("\t") for line in results[1:]]

    nucleotide_protein_seqs_df_domains = pd.DataFrame(data, columns=columns) # create a dataframe from the data
    nucleotide_protein_seqs_df_domains.drop(columns=["From"], inplace=True) # remove the "From" column
    nucleotide_protein_seqs_df_domains.columns = ['primary_id', 'pdb', 'domain', 'gene3d', 'family'] # rename the columns

    nucleotide_protein_seqs_df = pd.merge(nucleotide_protein_seqs_df, nucleotide_protein_seqs_df_domains, on="primary_id")

    sequnceDMF = sequenceDomainFetcher() # initilize the sequence domain fetcher
    # apply the assign_foldclass method to the dataframe in parallel to get the domain (and fold class) information
    nucleotide_protein_seqs_df.parallel_apply(lambda data: sequnceDMF.assign_foldclass(data=data, output_path=args.output_path, name=args.name), axis=1)

if __name__ == "__main__":
    main()