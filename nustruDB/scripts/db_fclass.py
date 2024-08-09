#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import shutil
import argparse

import logging
import pandas as pd

from pathlib import Path

from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True, nb_workers=10, verbose=0) 

from fetchUPTSV import UniProtFetcher

scop_fold_classes = {"1000000":"All alpha proteins", "1000001":"All beta proteins", "1000002":"Alpha and beta proteins (a/b)", 
                "1000003":"Alpha and beta proteins (a+b)", "1000004":"Small proteins"}

cath_fold_classes = {"1":"Mainly Alpha", "2":"Mainly Beta", "3":"Alpha Beta", "4":"Few Secondary Structures", "6":"Special"}

def extract_domains(text):
    # Regular expression to match the domain information
    domain_pattern = re.compile(r'DOMAIN (\d+)\.\.(\d+); /note="([^"]+)"; /evidence="([^"]+)"')
    
    # Find all matches in the text
    matches = domain_pattern.findall(text)
    
    # Create a list of lists for each domain
    domains = []
    for match in matches:
        start, end, note, evidence = match
        domains.append([start, end, note, evidence])
    
    return domains

def assign_foldclass(data, output_path, name, cfc= cath_fold_classes):
    # first try to assign domains that are annotated in the SCOP
    # hit = scopdb[scopdb["SF-UNIID"] == f"{data["primary_id"]}"]["SF-UNIREG"]

    # domains = data['domain'].split(' ')
    # domains = [list(domain) for k,domain in itertools.groupby(domains, lambda split: split=='DOMAIN') if not k]
    domains = extract_domains(data['domain'])
    
    gene3d = data['gene3d'].replace('"','').split(';')
    gene3d = [gene3d[i:i + 3] for i in range(0, len(gene3d), 3)]
    
    # supfam = data['supfam'].replace('"','').split(';')
    # supfam = [supfam[i:i + 3] for i in range(0, len(supfam), 3)]
    
    if len(domains) > 0 and len(domains[0]) > 0: 
        
        cath_fclass = []
        # supfam_fclass = []
        
        if len(gene3d[0]) > 0:
            for gene in gene3d:
                if len(gene) > 1:
                    cath_fclass.append(cfc.get(gene[0][0], "Unknown"))
        
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
        description="Fetch the annotated domains and fold classed for a protein."
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
    nucleotide_protein_seqs_df = pd.read_csv(args.input_file, nrows=100000)

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


    nucleotide_protein_seqs_df.parallel_apply(lambda data: assign_foldclass(data=data, output_path=args.output_path, name=args.name), axis=1)

if __name__ == "__main__":
    main()