#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import logging
import pandas as pd

from pathlib import Path
from Bio.Seq import Seq

from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True, nb_workers=8)

def nucleotide_to_protein(data, output_path, name):
    if data["nucleotide_sequence"][0:3] == "ATG":
        if len(data["nucleotide_sequence"]) % 3 == 0:
            if data["nucleotide_sequence"].count('A') + data["nucleotide_sequence"].count('T') + data["nucleotide_sequence"].count('G') + data["nucleotide_sequence"].count('C') == len(data["nucleotide_sequence"]):
                translated_sequence = Seq(data["nucleotide_sequence"][0:-3]).translate()
                
                if translated_sequence == data["protein_sequence"]:
                    with open(f'{output_path}/{name}_protein.fasta', 'a') as f, open(f'{output_path}/{name}_nucleotide.fasta', 'a') as f2:
                        f.write(f">{data["primary_id"]}\n{data["protein_sequence"]}\n")
                        f2.write(f">{data["primary_id"]}\n{data["nucleotide_sequence"]}\n")
                        
                    new_data = pd.DataFrame([{'source': data['source'], 'primary_id': data['primary_id'], 'gene_name': data['gene_name'], 'organism': data['organism'], 'expression_system': data['expression_system'],
                            'protein_sequence': data['protein_sequence'], 'nucleotide_id': data['nucleotide_id'], 'nucleotide_sequence': data['nucleotide_sequence'], 'bfactor_or_plddt': data['bfactor_or_plddt'], 'secondary_structure': data['secondary_structure']}])
                    new_data.to_csv(f'{output_path}/{name}.csv', mode='a', index=False, header=False)
                    del data, new_data

        else:
            logging.error(f"Error: {data["primary_id"]} nucleotide sequence does not translate to protein sequence.")
            pass
    else:
        logging.error(f"Error: {data["nucleotide_sequence"][0:3]} nucleotide sequence does not start with start codon.")
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
        '-w', '--overwrite', action="store_true", dest="overwrite", required=False, default=False,
        help='If file name already exists, overwrite it. Default is False.' 
    )
    args = parser.parse_args()

    
    logging.basicConfig(filename=f'{args.output_path}/{args.name}.log',
            filemode='w',
            format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
            datefmt='%H:%M:%S',
            level=logging.ERROR)

    
    df = pd.read_csv(args.input_file)
    
    if not Path(f'{args.output_path}/{args.name}.csv').exists() or args.overwrite:
        with open(f'{args.output_path}/{args.name}.csv', mode='w') as f:
            if "secondary_structure" in df.columns:
                f.write('source,primary_id,gene_name,organism,expression_system,protein_sequence,nucleotide_id,nucleotide_sequence,bfactor_or_plddt,secondary_structure\n')
            else:
                f.write('source,primary_id,gene_name,organism,expression_system,protein_sequence,nucleotide_id,nucleotide_sequence\n')
    else:
        print(f"Error: {args.output_path} already exists. Use -w to overwrite.")
        exit(1)
        
    df.parallel_apply(lambda data: nucleotide_to_protein(data=data, output_path=args.output_path, name=args.name), axis=1)
    
if __name__ == '__main__':
    main()
