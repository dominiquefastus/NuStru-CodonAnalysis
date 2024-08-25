#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
INFO: Pandarallel will run on 10 workers.
INFO: Pandarallel will use standard multiprocessing data transfer (pipe) to transfer data between the main process and workers.

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
  
Script to retrieve nucleotide sequences from uniprot IDs in a csv/ tsv file. The script will create a log file and a new csv file.
Example: python fastUPmapNT.py -i Example/examples_nustruDB/example_uniprotList.tsv -o . -n example_uniprotList_nustru [-w]
"""

import argparse

import pandas as pd 
import numpy as np
pd.options.mode.chained_assignment = None
pd.set_option('future.no_silent_downcasting', True)

from Bio import Entrez
from Bio import SeqIO

from pathlib import Path
from typing import List
import logging

# set the number of workers for parallel processing to 10 (maximum for API requests limit per second)
# if the computer has less cores, all cores will be used
from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True, nb_workers=10, verbose=0)

# lists to store cds and protein ids for ncbi
overall_cds = []
overall_protein_ids = []

class UPmapperNT:
    def __init__(self, mail, api_key) -> None:
        self.mail = mail
        self.api_key = api_key
            
    def retrieve_nucleotide_seq(self, entryID=None, protein_id=None):
        """Retrieve nucleotide sequence from NCBI nucleotide database"""
        Entrez.email = self.mail
        Entrez.api_key = self.api_key
        
        try:
            # retrieve nucleotide sequence from NCBI nucleotide database (nuccore)
            # retmode="text" returns the sequence in plain text and rettype="fasta_cds_na" returns the coding sequence
            with Entrez.efetch(db="nuccore", rettype="fasta_cds_na", retmode="text", id=entryID) as handle:
                # loop over the records in the fasta file and check if the protein id is in the description
                # sometimes the protein id has a wrong numbering, so we check for the correct one as well
                for seq_record in SeqIO.parse(handle, "fasta"):
                    if protein_id in seq_record.description or protein_id.replace('.2', '.1') in seq_record.description:
                        # assign the header and sequence to variables
                        head = seq_record.id + seq_record.description
                        sequence = f'{entryID}:{seq_record.seq}'
        except:
            logging.error(f"Error: Nucleotide sequence {entryID} not found.")
            pass

        return head, sequence

    def retrieve_protein_seq(self, entryID=None):
        """Retrieve protein sequence from NCBI protein database"""
        Entrez.email = self.mail
        Entrez.api_key = self.api_key
        
        try:
            try:
                # retrieve protein sequence from NCBI protein database (protein)
                # retmode="text" returns the sequence in plain text
                with Entrez.efetch(db="protein", rettype="fasta", retmode="text", id=entryID) as handle:
                    for seq_record in SeqIO.parse(handle, "fasta"):
                        sequence = seq_record.seq
            except:
                # check if the protein is a .1 or .2 version and retrieve the other version
                # so we try both versions to fetch the protein sequence
                if ".2" in entryID:
                    entryID = entryID.replace('.2',".1")
                else:
                    entryID = entryID.replace('.1',".2")
                with Entrez.efetch(db="protein", rettype="fasta", retmode="text", id=entryID) as handle:
                    # loop over the records in the fasta file and check if the protein id is in the description
                    for seq_record in SeqIO.parse(handle, "fasta"):
                        sequence = seq_record.seq
        except:
            logging.error(f"Error: Protein sequence {entryID} not found.")
            pass
                    
        return sequence

    def get_cds(self, data, output_path, name):
        """Get CDS from uniprot IDs"""
        try:
            # split the nucleotide id string in the dataframe and remove parantheses and whitespaces
            id = str(data['nucleotide_id']).strip().replace('"', '').replace(' ', '') # plain id
            id_list = id.split(';') # split the id string into a list by the semicolon
            del id_list[-1] # remove the last element of the list, which is an empty string
            
            # get the cds ids, protein ids and status from the id list by iterating over the list with different step sizes
            cds_ids = id_list[0::4]
            protein_ids = id_list[1::4]
            status = id_list[2::4]
            
            try:
                # loop over the cds ids, protein ids and status to get the correct cds and protein id
                for cds_id, protein_id, status in zip(cds_ids, protein_ids, status):
                    try:
                        # only if the status is empty, meaning no conflicts in the mapping, the cds and protein id are correct
                        if status == '-': 
                            # check if the protein sequence from uniprot matches the annotated protein sequence from ncbi           
                            if self.retrieve_protein_seq(protein_id) == data['protein_sequence']:
                                # if the proteins match, retrieve the referring nucleotide coding sequence from ncbi
                                matched_id, matched_seq = self.retrieve_nucleotide_seq(cds_id, protein_id)
                                break  
                            else:
                                continue
                        else:
                            continue
                    except:
                        # errors can be frameshifts, sequence conflicts, wrong translation, wrong cds id, etc.
                        logging.error(f"Error: {cds_id} and {protein_id} for {data['primary_id']} not found.")
                        continue
                
                # the nucleotide sequences can be created after filering the data             
                # with open(f'{output_path}/{name}.fasta', 'a') as f: # write the nucleotide sequence to a fasta file
                    # f.write(f">{data['primary_id']}|{matched_id} {data['organism']}\n{matched_seq.split(':')[1]}\n")

                data['nucleotide_sequence'] = matched_seq # assign the nucleotide sequence to the dataframe
                # build a new dataframe with the data and write it to a csv file (append mode)
                # avoiding dataloss during parallel processing (concurrent writes to the same file)
                data[['nucleotide_id','nucleotide_sequence']] = data['nucleotide_sequence'].split(':')
                new_data = pd.DataFrame({'source': 'uniprot', 'primary_id': data['primary_id'], 'gene_name': data['gene_name'], 'organism': data['organism'], 'expression_system': data['expression_system'],
                                        'protein_sequence': data['protein_sequence'], 'nucleotide_id': data['nucleotide_id'], 'nucleotide_sequence': data['nucleotide_sequence']}, index=[0])
                new_data.to_csv(f'{output_path}/{name}.csv', mode='a', index=False, header=False)
                logging.info(f"Success: {data['primary_id']} nucleotide sequence retrieved.")
                
                # delete the dataframes to free memory
                del data
                del new_data
                
            except:
                logging.error(f"Error: {cds_id} and {protein_id} for {data['primary_id']} not found.")
                pass
            
        except:
            logging.error(f"Error: CDS for {data['primary_id']} not available.")
            pass

def main():
    parser = argparse.ArgumentParser(
        prog='fastUPmapNT.py',
        description="Retrieve nucleotide sequences from uniprot IDs in a csv file."
    )
    parser.add_argument(
        '-i', '--input', type=str, dest="input_file", required=True,
        help='Input file with uniprot IDs in a csv format.'
    )
    parser.add_argument( 
        '-o', '--output', type=str, dest="output_path", required=True,
        help='Output path or directory to store the log file and the data.'
    )
    parser.add_argument(
        '-n', '--name', type=str, dest="name", required=True,
        help='Name of the output files and log file.'
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
    args = parser.parse_args()
    
    # create a log file
    logging.basicConfig(filename=f'{args.output_path}/{args.name}.log',
                    filemode='a',
                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.DEBUG)
    
    # read the input file and replace empty strings with NaN
    nucleotide_protein_seqs_df = pd.read_csv(args.input_file, sep='\t', header=0, dtype={"Subcellular location [CC]": object, "Alternative sequence": object})
    nucleotide_protein_seqs_df = nucleotide_protein_seqs_df.replace(r'^\s*$', np.nan, regex=True)

    nucleotide_protein_seqs_df = nucleotide_protein_seqs_df[['Entry', 'Gene Names (primary)', 'Organism', 'Subcellular location [CC]', 'Sequence', 'EMBL']]
    nucleotide_protein_seqs_df.columns = ['primary_id', 'gene_name', 'organism', 'expression_system', 'protein_sequence', 'nucleotide_id']
    
    # create a csv file with the header (columns)
    if not Path(f'{args.output_path}/{args.name}.csv').exists() or args.overwrite:
        with open(f'{args.output_path}/{args.name}.csv', mode='w') as f:
            f.write('source,primary_id,gene_name,organism,expression_system,protein_sequence,nucleotide_id,nucleotide_sequence\n')
    else:
        logging.error(f"Error: {args.output_path}/{args.name}.csv already exists.")
        exit(1)
    
    # create an instance of the UPmapperNT class
    mapper = UPmapperNT(mail=args.api_mail, api_key=args.api_key)
    
    # apply the nucleotide sequence retrieval function in parallel to the dataframe using pandarallel
    nucleotide_protein_seqs_df.parallel_apply(lambda data: mapper.get_cds(data=data, output_path=args.output_path, name=args.name), axis=1)
    
if __name__ == '__main__':
    main()