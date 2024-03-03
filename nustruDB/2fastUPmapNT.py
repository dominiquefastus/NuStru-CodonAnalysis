#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse

import pandas as pd 
import numpy as np
pd.options.mode.chained_assignment = None
pd.set_option('future.no_silent_downcasting', True)

from Bio import Entrez
from Bio import SeqIO

from typing import List
import logging

from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True)

overall_cds = []
overall_protein_ids = []

def retrieve_nucleotide_seq(entryID=None, protein_id=None):
    Entrez.email = "dominique_philip.fastus.5766@student.lu.se"    
    

    with Entrez.efetch(db="nuccore", rettype="fasta_cds_na", retmode="text", id=entryID) as handle:
        for seq_record in SeqIO.parse(handle, "fasta"):
            if protein_id in seq_record.description or protein_id.replace('.2', '.1') in seq_record.description:
                head = seq_record.id + seq_record.description
                sequence = f'{entryID}:{seq_record.seq}'

    return head, sequence

def retrieve_protein_seq(entryID=None):
    Entrez.email = "dominique_philip.fastus.5766@student.lu.se"    
    
    try:
        with Entrez.efetch(db="protein", rettype="fasta", retmode="text", id=entryID) as handle:
            for seq_record in SeqIO.parse(handle, "fasta"):
                sequence = seq_record.seq
    except:
        if ".2" in entryID:
            entryID = entryID.replace('.2',".1")
        else:
            entryID = entryID.replace('.1',".2")
        with Entrez.efetch(db="protein", rettype="fasta", retmode="text", id=entryID) as handle:
            for seq_record in SeqIO.parse(handle, "fasta"):
                sequence = seq_record.seq
                
    return sequence

def get_cds(data, output_path, name):
    id = str(data['nucleotide_id']).strip().replace('"', '').replace(' ', '')
    id_list = id.split(';')
    del id_list[-1]
    
    cds_ids = id_list[0::4]
    protein_ids = id_list[1::4]
    status = id_list[2::4]
    
    try:
        
        for cds_id, protein_id, status in zip(cds_ids, protein_ids, status):
            try:
                if status == '-':            
                    if retrieve_protein_seq(protein_id) == data['protein_sequence']:
                        
                        matched_id, matched_seq = retrieve_nucleotide_seq(cds_id, protein_id)
                        break  
                    else:
                        continue
                else:
                    continue
            except:
                logging.error(f"Error: {cds_id} and {protein_id} for {data['primary_id']} not found.")
                continue
                        
        with open(f'{output_path}/{name}.fasta', 'a') as f:
            f.write(f'>{matched_id} {data['organism']}\n{matched_seq.split(':')[1]}\n')

        data['nucleotide_sequence'] = matched_seq
        data[['nucleotide_id','nucleotide_sequence']] = data['nucleotide_sequence'].split(':')
        new_data = pd.DataFrame({'primary_id': data['primary_id'], 'gene_name': data['gene_name'], 'organism': data['organism'], 'expression_system': data['expression_system'],
                                'protein_sequence': data['protein_sequence'], 'nucleotide_id': data['nucleotide_id'], 'nucleotide_sequence': data['nucleotide_sequence']}, index=[0])
        new_data.to_csv(f'{output_path}/{name}.csv', mode='a', index=False, header=False)
        del data
        del new_data
        
    except:
        logging.error(f"Error: {cds_id} and {protein_id} for {data['primary_id']} not found.")
        pass

def main():
    
    parser = argparse.ArgumentParser(
        prog='2fastUPmapNT.py',
        description="Retrieve nucleotide sequences from uniprot IDs."
    )
    parser.add_argument(
        '-i', '--input', type=str, dest="input_file", required=True,
        help='Input file with uniprot IDs.'
    )
    parser.add_argument( 
        '-o', '--output', type=str, dest="output_path", required=True,
        help='Output file with nucleotide sequences.'
    )
    parser.add_argument(
        '-n', '--name', type=str, dest="name", required=True,
        help='Name of the output files and log file.'
    )
    args = parser.parse_args()
    
    logging.basicConfig(filename=f'{args.output_path}/{args.name}.log',
                    filemode='a',
                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.DEBUG)
    
    df = pd.read_csv(args.input_file, sep='\t', header=0, nrows=20, dtype={"Subcellular location [CC]": object, "Alternative sequence": object})
    df = df.replace(r'^\s*$', np.nan, regex=True)

    nucleotide_protein_seqs_df = df[['Entry', 'Gene Names (primary)', 'Organism', 'Subcellular location [CC]', 'Sequence', 'EMBL']]
    nucleotide_protein_seqs_df.columns = ['primary_id', 'gene_name', 'organism', 'expression_system', 'protein_sequence', 'nucleotide_id']
    
    with open(f'{args.output_path}/{args.name}.csv', mode='w') as f:
        f.write('primary_id,gene_name,organism,expression_system,protein_sequence,nucleotide_id,nucleotide_sequence\n')
        
    nucleotide_protein_seqs_df.parallel_apply(lambda data: get_cds(data=data, output_path=args.output_path, name=args.name), axis=1)
    
if __name__ == '__main__':
    main()
