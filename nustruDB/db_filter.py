#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import logging
import pandas as pd

from Bio.Seq import Seq

def nucleotide_to_protein(columns, row, output_path, name):

    if row.nucleotide_sequence[0:3] == "ATG":
        if len(row.nucleotide_sequence) % 3 == 0:
            if row.nucleotide_sequence.count('A') + row.nucleotide_sequence.count('T') + row.nucleotide_sequence.count('G') + row.nucleotide_sequence.count('C') == len(row.nucleotide_sequence):
                translated_sequence = (Seq(row.nucleotide_sequence[0:-3]).translate())
                
                if translated_sequence == row.protein_sequence:
                    with open(f'{output_path}/{name}_protein.fasta', 'a') as f, open(f'{output_path}/{name}_nucleotide.fasta', 'a') as f2:
                        f.write(f">{row.primary_id}\n{row.protein_sequence}\n")
                        f2.write(f">{row.primary_id}\n{row.nucleotide_sequence}\n")
                        
                    new_data = pd.DataFrame(dict(zip(columns, row.tolist())), index=[0])
                    new_data.to_csv(f'{output_path}/{name}.csv', mode='a', index=False, header=False)

        else:
            logging.error(f"Error: {row.primary_id} nucleotide sequence does not translate to protein sequence.")
            pass
    else:
        logging.error(f"Error: {row.nucleotide_sequence[0:3]} nucleotide sequence does not start with start codon.")
        pass
        
def main():

    parser = argparse.ArgumentParser(
        prog='db_filter.py',
        description="Filter the csv databases to reduce rendundancy, check cds and right translation and further improve the data contained."
    )
    parser.add_argument(
        '-i', '--input', type=str, dest="input_file", required=True,
        help='Input file of csv formatted uniprot and/ or pdb entries.'
    )
    parser.add_argument( 
        '-o', '--output', type=str, dest="output_path", required=True,
        help='Output to store the new filtered and reduced csv.'
    )
    parser.add_argument(
        '-n', '--name', type=str, dest="name", required=True,
        help='Name of the output csv file and log file.'
    )
    args = parser.parse_args()

    
    logging.basicConfig(filename=f'{args.output_path}/{args.name}.log',
            filemode='w',
            format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
            datefmt='%H:%M:%S',
            level=logging.ERROR)

    
    df = pd.read_csv('/Users/dominiquefastus/master_project/NuStru/nustruDB/NEW_ECOLI_FULL_uniprot_02_sec_struc.csv')
    
    with open(f'{args.output_path}/{args.name}.csv', mode='w') as f:
        f.write('source,primary_id,gene_name,organism,expression_srow.protein_sequenceucleotirow.nucleotide_sequence')

    for row in df.iterrows():
        nucleotide_to_protein(columns=df.columns, row=row[1], output_path=args.output_path, name=args.name)
    
if __name__ == '__main__':
    main()
