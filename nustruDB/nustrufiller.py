#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import argparse
import asyncio
import logging
from pathlib import Path
import shutil

import pandas as pd 
import numpy as np
pd.options.mode.chained_assignment = None
pd.set_option('future.no_silent_downcasting', True) 

from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True, nb_workers=10, verbose=0)

from scripts.fetchINPRO import InterProFetcher
from scripts.fetchUPTSV import UniProtFetcher
from scripts.fastUPmapNT import UPmapperNT
from scripts.db_filter import filter_sequences
from scripts.db_fetch import sequenceFeautureFetcher

async def main():
    parser = argparse.ArgumentParser(
            prog='nustrufiller.py',
            description="Complete pipeline to fetch the data for uniprot IDs or interpro IDs (complete protein families)."
        )
    parser.add_argument(
        '-i', '--input', type=str, dest="input", required=True,
        help='Single uniprotID or interproID, or file with comma seperated ids.'
        )
    parser.add_argument( 
        '-o', '--output', type=str, dest="output_path", required=True,
        help='Output path or directory to store the log file and the data.'
        )
    parser.add_argument(
        '-n', '--name', type=str, dest="name", required=False,
        help='Name of the output files and log file.'
        )
    parser.add_argument(
        '-u', '--unique', type=str, dest="unique", required=False,
        help='pssibility to drop duplicate. To keep duplicates use None. Default: organisms' 
    )
    parser.add_argument(
        '-w', '--overwrite', action="store_true", dest="overwrite", required=False, default=False,
        help='If file name already exists, overwrite it. Default is False.' 
        )
    parser.add_argument(
        '-m', '--mail', type=str, dest="api_mail", required=True, default=False,
        help='Provide mail for ncbi api account: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/.' 
    )
    parser.add_argument(
        '-k', '--key', type=str, dest="api_key", required=True, default=False,
        help='Provide api key from ncbi api account: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/.' 
    )
    parser.add_argument(
        '-d', '--download', action="store_true", dest="download", required=False, default=False,
        help='Download the structure files.'
    )
    args = parser.parse_args()
    
    if args.name is None:
        args.name = Path(args.input).stem if Path(args.input).exists() else args.input
    
    dir_name = Path(args.input).stem if Path(args.input).exists() else args.input
    Path(f'{args.output_path}/{dir_name}_nustruDB').mkdir(parents=True, exist_ok=True)
        
    logging.basicConfig(filename=f'{args.output_path}/{dir_name}_nustruDB/{args.name}.log',
                    filemode='a',
                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.ERROR)
    
    ##############################################################################################################
    print("--------------------------------------------")
    if Path(args.input).exists():
        dir_name = Path(args.input).stem
        with open(args.input, 'r') as file:
            ids = file.read().split(',')
        
        if ids[0].startswith('IPR'):
            print(f"Fetching the uniprot IDs for {len(ids)} protein families...")
            
            for id in ids:
                IPfetcher = InterProFetcher()
                await IPfetcher.run(family_id=id, output_path=f"{args.output_path}/{dir_name}_nustruDB")
                
                with open(f"{args.output_path}/{args.family_id}_nustruDB/{args.family_id}_fetched.txt", "r") as file:
                    uniprotIDs = file.read().split(',')
        
        else:
            with open(args.input, 'r') as file:
                uniprotIDs = file.read().splitlines()
            
    else:
        dir_name = args.input
        if args.input.startswith('IPR'):
            print(f"Fetching the uniprot IDs for protein family {args.input}...")
            IPfetcher = InterProFetcher()
            await IPfetcher.run(family_id=args.input, output_path=f"{args.output_path}/{dir_name}_nustruDB")
            
            with open(f"{args.output_path}/{dir_name}_nustruDB/{args.input}_fetched.txt", "r") as file:
                uniprotIDs = file.read().split(',')
        
        else:
            uniprotIDs = args.input.split(' ')
     
    ##############################################################################################################  
    print("--------------------------------------------")     
    print(f"Fetching the data for {len(uniprotIDs)} uniprot IDs...")
    
    UPfetcher = UniProtFetcher()
        # submit the ID mapping job and get the job ID
    job_id = UPfetcher.submit_id_mapping(from_db="UniProtKB_AC-ID", to_db="UniProtKB", ids=uniprotIDs)
    
    # check if the ID mapping results are ready and get the results
    if UPfetcher.check_id_mapping_results_ready(job_id):
        link = UPfetcher.get_id_mapping_results_link(job_id)
        # get the ID mapping results from the search URL and store them in a tsv file with the defined features as columns
        results = UPfetcher.get_id_mapping_results_search(link + "?format=tsv" + "&fields=accession%2Cgene_primary%2Corganism_name%2Ccc_subcellular_location%2Csequence%2Cxref_embl_full%2Cxref_ccds%2Cxref_refseq_full%2Ccc_alternative_products%2Cft_var_seq") 

    # store the results in a tsv file with the defined features as columns 
    columns = results[0].split("\t") 
    data = [line.split("\t") for line in results[1:]]
    
    df = pd.DataFrame(data, columns=columns) # create a dataframe from the data
    df.drop(columns=["From"], inplace=True) # remove the "From" column
    
    # check if the output path exists and write the dataframe to a tsv file
    if Path(args.output_path).exists() or args.overwrite:
        df.to_csv(f"{args.output_path}/{dir_name}_nustruDB/{args.name}.tsv", sep="\t", index=False)
    else:
        print(f"Error: {args.output_path}/{dir_name}_nustruDB/ already exists. Use -w to overwrite.")
        exit(1)
        
    ##############################################################################################################
    # get the nucleotide sequences from the uniprot IDs
    print("--------------------------------------------")
    print(f"Fetching the nucleotide sequences...")
    df = pd.read_csv(f"{args.output_path}/{dir_name}_nustruDB/{args.name}.tsv", sep='\t', header=0, dtype={"Subcellular location [CC]": object, "Alternative sequence": object})
    df = df.replace(r'^\s*$', np.nan, regex=True)

    nucleotide_protein_seqs_df = df[['Entry', 'Gene Names (primary)', 'Organism', 'Subcellular location [CC]', 'Sequence', 'EMBL']]
    nucleotide_protein_seqs_df.columns = ['primary_id', 'gene_name', 'organism', 'expression_system', 'protein_sequence', 'nucleotide_id']

    if not Path(f'{args.output_path}/{dir_name}_nustruDB/{args.name}.csv').exists() or args.overwrite:
        with open(f'{args.output_path}/{dir_name}_nustruDB/{args.name}.csv', mode='w') as f:
            f.write('source,primary_id,gene_name,organism,expression_system,protein_sequence,nucleotide_id,nucleotide_sequence\n')
    else:
        logging.error(f"Error: {args.output_path}/{dir_name}_nustruDB/{args.name}.csv already exists.")
        exit(1)
        
    mapper = UPmapperNT(mail=args.api_mail, api_key=args.api_key)
    nucleotide_protein_seqs_df.parallel_apply(lambda data: mapper.get_cds(data=data, output_path=f"{args.output_path}/{dir_name}_nustruDB", name=args.name), axis=1)
    
    ##############################################################################################################
    # filter the data in parallel and remove duplicates by a column from the new csv file
    print("--------------------------------------------")
    print(f"Filtering the data...")
    # read the csv file as input
    df = pd.read_csv(f"{args.output_path}/{dir_name}_nustruDB/{args.name}.csv", index_col=False)
    
    # create the new csv file with the respective columns if it does not exist or the overwrite flag is set
    # write the header to the new csv file depending on the columns in the data
    if not Path(f"{args.output_path}/{dir_name}_nustruDB/{args.name}_filtered.csv").exists() or args.overwrite:
        with open(f"{args.output_path}/{dir_name}_nustruDB/{args.name}_filtered.csv", mode='w') as f:
            if "secondary_structure" in df.columns:
                f.write('source,primary_id,gene_name,organism,expression_system,protein_sequence,nucleotide_id,nucleotide_sequence,bfactor_or_plddt,secondary_structure\n')
            else:
                f.write('source,primary_id,gene_name,organism,expression_system,protein_sequence,nucleotide_id,nucleotide_sequence\n')
    else:
        print(f"Error: {args.output_path}/{dir_name}_nustruDB already exists. Use -w to overwrite.")
        exit(1)
        
    # apply the nucleotide_to_protein function to filter the data in parallel
    df.parallel_apply(lambda data: filter_sequences(data=data, output_path=f"{args.output_path}/{dir_name}_nustruDB/", 
                                                    name=f"{args.name}_filtered", prot_fasta=False, nuc_fasta=False), axis=1)
    
    # remove duplicates by a column from the new csv file and overwrite the new csv file if the flag is set
    if args.unique is not None:
        nustrudb = pd.read_csv(f"{args.output_path}/{dir_name}_nustruDB/{args.name}_filtered.csv")
        unique_df = nustrudb.drop_duplicates(subset=args.unique)
        unique_df.to_csv(f"{args.output_path}/{dir_name}_nustruDB/{args.name}_filtered.csv", index=False, header=True)
        
    ##############################################################################################################
    print("--------------------------------------------")
    print(f"Fetching the structure files...")
    
    # read the input file
    nucleotide_protein_seqs_df = pd.read_csv(f"{args.output_path}/{dir_name}_nustruDB/{args.name}_filtered.csv", index_col=False)
    
    # create the new csv file with the respective columns if it does not exist or the overwrite flag is set
    if not Path(f"{args.output_path}/{dir_name}_nustruDB/{args.name}_filtered_secstru.csv").exists() or args.overwrite:
        with open(f"{args.output_path}/{dir_name}_nustruDB/{args.name}_filtered_secstru.csv", mode='w') as f:
            f.write('source,primary_id,gene_name,organism,expression_system,protein_sequence,nucleotide_id,nucleotide_sequence,bfactor_or_plddt,secondary_structure\n')
    else:
        print(f"Error: {args.output_path}/{dir_name}_nustruDB already exists. Use -w to overwrite.")
        exit(1)
    
    sequenceFF = sequenceFeautureFetcher() # initialize the sequence feature fetcher
    # apply the fetch_pdb_and_plddt function to get the features to the dataframe in parallel
    nucleotide_protein_seqs_df.parallel_apply(lambda data: sequenceFF.fetch_pdb_and_plddt(data=data, output_path=f"{args.output_path}/{dir_name}_nustruDB/", name=f"{args.name}_filtered_secstru", download=args.download), axis=1)
    
    # delete the cif and pdb directories if the download flag is set to False
    if not args.download:
        shutil.rmtree(f'{args.output_path}/{dir_name}_nustruDB/cif_files/')
        shutil.rmtree(f'{args.output_path}/{dir_name}_nustruDB/pdb_files/')
    
if __name__ == '__main__':
    asyncio.run(main())