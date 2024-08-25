#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
INFO: Pandarallel will run on 8 workers.
INFO: Pandarallel will use standard multiprocessing data transfer (pipe) to transfer data between the main process and workers.
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
  
Fetching the with the previous data fetching scripts can still contain wrong assigned sequences or translation mistakes. To address this
issue, the script will check if the nucleotide sequence can be translated to the related protein sequence. It also removes duplicates by a column.
The script will write the new csv file and log file in the output path. 
Example: python db_filter.py -i Example/examples_nustruDB/example_db_filtered.csv -o . -n nustruDB -u 
"""

import argparse
import logging
import pandas as pd

from pathlib import Path
from Bio.Seq import Seq

from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True, nb_workers=8, verbose=0)

def filter_sequences(data, output_path, name, create_fasta='none'):
    """Checks and writed matching nucleotide and protein sequence"""
    # check if the nucleotide sequence starts with start codon (normally ATG)
    if data['nucleotide_sequence'][0:3] == "ATG":
  
        # check if the nucleotide sequence consists of an even number of codons
        if len(data['nucleotide_sequence']) % 3 == 0:
            # check if the nucleotide sequence consists of only A, T, G, C
            if data['nucleotide_sequence'].count('A') + data['nucleotide_sequence'].count('T') + data['nucleotide_sequence'].count('G') + data['nucleotide_sequence'].count('C') == len(data['nucleotide_sequence']):
                # translate the nucleotide sequence to protein sequence without the stop codon
                
                if data['source'] == 'pdb':
                    # depending on the source, the nucleotide sequence is needed or not
                    translated_sequence = Seq(data['nucleotide_sequence']).translate()
                else:
                    translated_sequence = Seq(data['nucleotide_sequence'][0:-3]).translate()
       
                # check if the translated protein sequence is equal to the protein sequence in the data
                # then write the protein and nucleotide sequence to a fasta file
                if translated_sequence == data['protein_sequence']:
                    # check if a protein and nucleotide fasta files should be created
                    if create_fasta == 'protein' or create_fasta == 'all':
                        with open(f'{output_path}/{name}_protein.fasta', 'a') as f:
                            f.write(f">{data['primary_id']}\n{data['protein_sequence']}\n")
                    if create_fasta == 'nucleotide' or create_fasta == 'all':
                        with open(f'{output_path}/{name}_nucleotide.fasta', 'a') as f2:
                            f2.write(f">{data['primary_id']}\n{data['nucleotide_sequence']}\n")
                    
                    # write the data to a new csv file filtered entries and delete the data to free memory
                    # if the data contains secondary structure, write it to the new csv file    
                    print("Ture")
                    if "secondary_structure" in data.keys():
                        new_data = pd.DataFrame([{'source': data['source'], 'primary_id': data['primary_id'], 'gene_name': data['gene_name'], 'organism': data['organism'], 'expression_system': data['expression_system'],
                            'protein_sequence': data['protein_sequence'], 'nucleotide_id': data['nucleotide_id'], 'nucleotide_sequence': data['nucleotide_sequence'], 'bfactor_or_plddt': data['bfactor_or_plddt'], 'secondary_structure': data['secondary_structure']}])
                    else:
                        new_data = pd.DataFrame([{'source': data['source'], 'primary_id': data['primary_id'], 'gene_name': data['gene_name'], 'organism': data['organism'], 'expression_system': data['expression_system'],
                                'protein_sequence': data['protein_sequence'], 'nucleotide_id': data['nucleotide_id'], 'nucleotide_sequence': data['nucleotide_sequence']}])
                    new_data.to_csv(f'{output_path}/{name}.csv', mode='a', index=False, header=False)
                    del data, new_data

        else:
            logging.error(f"Error: {data['primary_id']} nucleotide sequence does not translate to protein sequence.")
            pass
    else:
        logging.error(f"Error: {data['nucleotide_sequence'][0:3]} nucleotide sequence does not start with start codon.")
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
    parser.add_argument(
        '-u', '--unique', type=str, dest="unique", required=False,
        help='pssibility to drop duplicate. To keep duplicates use None. Default: organisms' 
    ) 
    parser.add_argument(
        '--create-fasta', choices=['protein', 'nucleotide', 'all', 'none'], dest="create_fasta", default='none',
        help="Create a fasta file with the nucleotide sequences, protein sequences or both."
    )
    parser.add_argument(
        '-w', '--overwrite', action="store_true", dest="overwrite", required=False, default=False,
        help='If file name already exists, overwrite it. Default is False.' 
    )
    args = parser.parse_args()

    logging.basicConfig(filename=f'{args.output_path}/{args.name}.log',
            filemode='w',
            format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
            datefmt='%H:%M:%S',
            level=logging.ERROR)

    # read the csv file as input
    nucleotide_protein_seqs_df = pd.read_csv(args.input_file)
    
    # create the new csv file with the respective columns if it does not exist or the overwrite flag is set
    # write the header to the new csv file depending on the columns in the data
    if not Path(f'{args.output_path}/{args.name}.csv').exists() or args.overwrite:
        with open(f'{args.output_path}/{args.name}.csv', mode='w') as f:
            if "secondary_structure" in nucleotide_protein_seqs_df.columns:
                f.write('source,primary_id,gene_name,organism,expression_system,protein_sequence,nucleotide_id,nucleotide_sequence,bfactor_or_plddt,secondary_structure\n')
            else:
                f.write('source,primary_id,gene_name,organism,expression_system,protein_sequence,nucleotide_id,nucleotide_sequence\n')
    else:
        print(f"Error: {args.output_path} already exists. Use -w to overwrite.")
        exit(1)
        
    # apply the nucleotide_to_protein function to filter the data in parallel
    nucleotide_protein_seqs_df.parallel_apply(lambda data: filter_sequences(data=data, output_path=args.output_path, name=args.name), axis=1)
    
    # remove duplicates by a column from the new csv file and overwrite the new csv file if the flag is set
    if args.unique is not None:
        nustrudb = pd.read_csv(f'{args.output_path}/{args.name}.csv')
        unique_df = nustrudb.drop_duplicates(subset=args.unique)
        unique_df.to_csv(f'{args.output_path}/{args.name}.csv', index=False, header=True)
    
if __name__ == '__main__':
    main()
