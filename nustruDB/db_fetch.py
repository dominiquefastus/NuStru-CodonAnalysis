#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import numpy as np

from pathlib import Path

from urllib.request import urlretrieve
from urllib.error import HTTPError

from biopandas.pdb import PandasPdb 

from Bio.PDB import PDBParser
from Bio.PDB import MMCIFParser
from Bio.PDB.DSSP import DSSP

from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True)

def fetch_pdb_and_plddt(data, output_path, name):
    pdbp = PDBParser(QUIET=True)
    cifp = MMCIFParser(QUIET=True)
    
    Path(f'{output_path}/pdb_files').mkdir(parents=True, exist_ok=True)
    Path(f'{output_path}/cif_files').mkdir(parents=True, exist_ok=True)
    
    if data['source'] == 'uniprot':
        print(data['primary_id'])
        ppdb = PandasPdb().fetch_pdb(uniprot_id=data['primary_id'], source="alphafold2-v4")
        
        af2_version = 4
        urlretrieve(f"https://alphafold.ebi.ac.uk/files/AF-{data['primary_id']}-F1-model_v{af2_version}.cif", f"{output_path}/cif_files/{data['primary_id']}.cif")
        
        residue_plddt = ppdb.df["ATOM"][["b_factor","residue_number"]]


        residue_plddt = residue_plddt.groupby(["residue_number"]).mean().round(4)

        residue_plddt_dict = residue_plddt.to_dict()
        data['bfactor_or_plddt'] = residue_plddt_dict
        
    else:
        print(data['primary_id'])
        ppdb = PandasPdb().fetch_pdb(pdb_code=data['primary_id'], source="pdb")
        
        Path(f'{output_path}/pdb_files').mkdir(parents=True, exist_ok=True)
        urlretrieve(f"https://files.rcsb.org/download/{data["primary_id"]}.cif", f"{output_path}/cif_files/{data['primary_id']}.cif")
        
        residue_b_factor = ppdb.df["ATOM"][["b_factor","residue_number"]]
        
        residue_b_factor = residue_b_factor.groupby(["residue_number"]).mean().round(4)

        residue_b_factor_dict = residue_b_factor.to_dict()
        
        data['bfactor_or_plddt'] = residue_b_factor_dict
    

    ppdb.to_pdb(path=f'{output_path}/pdb_files/{data['primary_id']}.pdb')

    structure = cifp.get_structure(data["primary_id"], f"{output_path}/cif_files/{data['primary_id']}.cif")
    model = structure[0]
    dssp = DSSP(model, f"{output_path}/{data['primary_id']}.cif")
    secondary_structure = "".join([dssp_data[2] for dssp_data in dssp])
    data['secondary_structure'] = secondary_structure

    new_data = pd.DataFrame({'source': data['source'], 'primary_id': data['primary_id'], 'gene_name': data['gene_name'], 'organism': data['organism'], 'expression_system': data['expression_system'],
                            'protein_sequence': data['protein_sequence'], 'nucleotide_id': data['nucleotide_id'], 'nucleotide_sequence': data['nucleotide_sequence'], 'bfactor_or_plddt': data['bfactor_or_plddt'], 'secondary_structure': data['secondary_structure']})
    new_data.to_csv(f'{output_path}/{name}.csv', mode='a', index=False, header=False)
    del new_data

def main():
    
    parser = argparse.ArgumentParser(
        prog='db_fetch.py',
        description="Fetch the pdb and map the plddt or bfactor to the protein sequence."
    )
    parser.add_argument(
        '-i', '--input', type=str, dest="input_file", required=True,
        help='Input file of csv formatted uniprot or pdb entries.'
    )
    parser.add_argument( 
        '-o', '--output', type=str, dest="output_path", required=True,
        help='Output to store the new csv with secondary structure information.'
    )
    parser.add_argument(
        '-n', '--name', type=str, dest="name", required=True,
        help='Name of the output files and log file.'
    )

    args = parser.parse_args()
    
    nucleotide_protein_seqs_df = pd.read_csv(args.input_file, index_col=False)
    
    with open(f'{args.output_path}/{args.name}.csv', mode='w') as f:
        f.write('source,primary_id,gene_name,organism,expression_system,protein_sequence,nucleotide_id,nucleotide_sequence,bfactor_or_plddt,secondary_structure\n')
        
    nucleotide_protein_seqs_df.apply(lambda data: fetch_pdb_and_plddt(data=data, output_path=args.output_path, name=args.name), axis=1)
    
if __name__ == "__main__":
    main()