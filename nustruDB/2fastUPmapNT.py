#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd 
import numpy as np
pd.options.mode.chained_assignment = None
pd.set_option('future.no_silent_downcasting', True)

from Bio import Entrez
from Bio import SeqIO

from typing import List
import logging
logging.basicConfig(filename='/Users/dominiquefastus/Downloads/log.log',
                    filemode='a',
                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.DEBUG)

from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True)

overall_cds = []
overall_protein_ids = []

def retrieve_nucleotide_seq(entryID=None, protein_id=None):
    Entrez.email = "dominique.fastus@biochemistry.lu.se"    
    

    with Entrez.efetch(db="nuccore", rettype="fasta_cds_na", retmode="text", id=entryID) as handle:
        for seq_record in SeqIO.parse(handle, "fasta"):
            if protein_id in seq_record.description or protein_id.replace('.2', '.1') in seq_record.description:
                head = seq_record.id + seq_record.description
                sequence = f'{entryID}:{seq_record.seq}'

    return head, sequence

def retrieve_protein_seq(entryID=None):
    Entrez.email = "dominique.fastus@biochemistry.lu.se"    
    
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

def get_cds(data):
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
                        
        with open('/Users/dominiquefastus/Downloads/sequences.fasta', 'a') as f:
            f.write(f'>{matched_id}\n{matched_seq}\n')

        data['nucleotide_sequence'] = matched_seq
        data[['nucleotide_id','nucleotide_sequence']] = data['nucleotide_sequence'].split(':')
        new_data = pd.DataFrame({'primary_id': data['primary_id'], 'gene_name': data['gene_name'], 'organism': data['organism'], 'expression_system': data['expression_system'],
                                'protein_sequence': data['protein_sequence'], 'nucleotide_id': data['nucleotide_id'], 'nucleotide_sequence': data['nucleotide_sequence']}, index=[0])
        new_data.to_csv('/Users/dominiquefastus/Downloads/nucleotide_protein_seqs.csv', mode='a', index=False, header=False)
        del data
        del new_data
        
    except:
        logging.error(f"Error: {cds_id} and {protein_id} for {data['primary_id']} not found.")
        pass

def main():
    
    df = pd.read_csv('/Users/dominiquefastus/master_project/NuStru/nustruDB/uniprotkb_taxonomy_id_562_2024_02_19.tsv', sep='\t', header=0, nrows=400, dtype={"Subcellular location [CC]": object,
                                                                                                                         "Alternative sequence": object})
    df = df.replace(r'^\s*$', np.nan, regex=True)

    nucleotide_protein_seqs_df = df[['Entry', 'Gene Names (primary)', 'Organism', 'Subcellular location [CC]', 'Sequence', 'EMBL']]
    nucleotide_protein_seqs_df.columns = ['primary_id', 'gene_name', 'organism', 'expression_system', 'protein_sequence', 'nucleotide_id']
    
    with open('/Users/dominiquefastus/Downloads/nucleotide_protein_seqs.csv', mode='w') as f:
        f.write('primary_id,gene_name,organism,expression_system,protein_sequence,nucleotide_id,nucleotide_sequence\n')
        
    nucleotide_protein_seqs_df.parallel_apply(get_cds, axis=1)
    
if __name__ == '__main__':
    main()
