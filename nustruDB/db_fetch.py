#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import shutil
import argparse

import logging
import pandas as pd

from pathlib import Path

import requests
from requests.adapters import HTTPAdapter, Retry

session = requests.Session()
retries = Retry(total=6, backoff_factor=0.2, status_forcelist=[ 502, 503, 504 ])
session.mount('https://', HTTPAdapter(max_retries=retries))

from biopandas.pdb import PandasPdb 

from Bio.PDB import MMCIFParser
cifp = MMCIFParser(QUIET=True)
    
from Bio.PDB.DSSP import DSSP
from Bio.SeqUtils import seq1


from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True)

def reassign_position(sequence, residue_position):
    substring = "".join(residue_position.values())
    res = re.search(substring, sequence)
    
    if res:
        position_range = range(res.start(), res.end())
        new_dict = {}
        number_dict = {}

        for position_seq, position_dict in zip(position_range, residue_position):
            new_dict[position_seq + 1] = residue_position[position_dict]
            number_dict[position_dict] = position_seq + 1

        return new_dict, number_dict
    else:
        print("Substring not found in the string.")
        return None

def fetch_features(data, name, chain_id, sequence, output_path, model_id):
    ppdb = PandasPdb().read_pdb(f"{output_path}/pdb_files/{model_id}.pdb")
     
    residue_position = ppdb.df["ATOM"][ppdb.df['ATOM']['chain_id'] == chain_id][["residue_number","residue_name"]]
    residue_position = residue_position.groupby(["residue_number"])["residue_name"].last().apply(seq1)
    residue_position_dict = residue_position.to_dict()
    residue_position_dict, number_dict = reassign_position(sequence=sequence, residue_position=residue_position_dict)
    
    residue_b_factor = ppdb.df["ATOM"][ppdb.df['ATOM']['chain_id'] == chain_id][["b_factor","residue_number"]]
    residue_b_factor = residue_b_factor.groupby(["residue_number"]).mean().round(4)['b_factor']
    residue_b_factor_dict = residue_b_factor.to_dict()
    
    updated_residue_b_factor_dict = {}
    for position in residue_b_factor_dict.keys():
        new_position = number_dict[position]
        updated_residue_b_factor_dict[new_position] = residue_b_factor[position]
    
    
    structure = cifp.get_structure(model_id, f"{output_path}/cif_files/{model_id}.cif")

    model = structure[0]
    dssp = DSSP(model, f"{output_path}/cif_files/{model_id}.cif")
    secondary_structure = [dssp_data[0:3:2] for dssp_data in dssp if dssp_data[0] in residue_position_dict.keys()]
    secondary_structure_dict = dict(secondary_structure)
    
    return updated_residue_b_factor_dict, secondary_structure_dict

def fetch_pdb_and_plddt(data, output_path, name, download=False):
    Path(f'{output_path}/pdb_files').mkdir(parents=True, exist_ok=True)
    Path(f'{output_path}/cif_files').mkdir(parents=True, exist_ok=True)

    try:
        if data['source'] == 'uniprot':
            model_id = data['primary_id']
            chain_id = 'A'
            
            try:
                af2_version = 4
                response_cif = session.get(f"https://alphafold.ebi.ac.uk/files/AF-{model_id}-F1-model_v{af2_version}.cif")
                response_pdb = session.get(f"https://alphafold.ebi.ac.uk/files/AF-{model_id}-F1-model_v{af2_version}.pdb")
                
                if "Error" not in response_cif.text and "Error" not in response_pdb.text:           
                    with open(f"{output_path}/cif_files/{model_id}.cif", "wb") as cif_file, open(f"{output_path}/pdb_files/{model_id}.pdb", "wb") as pdb_file:
                        cif_file.write(response_cif.content)
                        pdb_file.write(response_pdb.content)
    
                    residue_b_factor_dict,secondary_structure_dict = fetch_features(data=data, name=name, model_id=model_id, chain_id=chain_id, sequence=data['protein_sequence'], output_path=output_path)
                    
                else:
                    logging.error(f"Error: could not retrieve cif or file for {model_id}.")
                    pass
            except: 
                logging.error(f"Could not assign plddt to {model_id}.")
                pass
            
        else:
            model_id = data['primary_id'].replace('"','').split("_")[0]
            chain_id = data['primary_id'].replace('"','').split("_")[1]
            
            try:
                response_cif = session.get(f"https://files.rcsb.org/download/{model_id}.cif")
                response_pdb = session.get(f"https://files.rcsb.org/download/{model_id}.pdb")
                
                if "Error" not in response_cif.text and "Error" not in response_pdb.text:           
                    with open(f"{output_path}/cif_files/{model_id}.cif", "wb") as cif_file, open(f"{output_path}/pdb_files/{model_id}.pdb", "wb") as pdb_file:
                        cif_file.write(response_cif.content)
                        pdb_file.write(response_pdb.content)
                    
                    residue_b_factor_dict,secondary_structure_dict = fetch_features(data=data, name=name, model_id=model_id, chain_id=chain_id, sequence=data['protein_sequence'], output_path=output_path)
                    
                else:
                    logging.error(f"Error: could not retrieve cif or pdb file for {model_id}.")
                    pass
            except:
                logging.error(f"Could not assign plddt to {model_id}.")
                pass
        
        if not download:
            Path(f'{output_path}/pdb_files/{model_id}.pdb').unlink()
            
        data['bfactor_or_plddt'] = residue_b_factor_dict
        data['secondary_structure'] = secondary_structure_dict

        new_data = pd.DataFrame([{'source': data['source'], 'primary_id': data['primary_id'], 'gene_name': data['gene_name'], 'organism': data['organism'], 'expression_system': data['expression_system'],
                'protein_sequence': data['protein_sequence'], 'nucleotide_id': data['nucleotide_id'], 'nucleotide_sequence': data['nucleotide_sequence'], 'bfactor_or_plddt': data['bfactor_or_plddt'], 'secondary_structure': data['secondary_structure']}])
        new_data.to_csv(f'{output_path}/{name}.csv', mode='a', index=False, header=False)
        del data
        del new_data
        if not download:
            Path(f'{output_path}/cif_files/{model_id}.cif').unlink()
            
    except:
        logging.error(f"Error: {data['primary_id']} not found.")
        pass

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
    parser.add_argument(
        '-d', '--download', action="store_true", dest="download", required=False, default=False,
        help='Download the pdb files.'
    )
    args = parser.parse_args()
    
    logging.basicConfig(filename=f'{args.output_path}/{args.name}.log',
            filemode='a',
            format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
            datefmt='%H:%M:%S',
            level=logging.ERROR)
    
    nucleotide_protein_seqs_df = pd.read_csv(args.input_file, index_col=False)
    
    with open(f'{args.output_path}/{args.name}.csv', mode='w') as f:
        f.write('source,primary_id,gene_name,organism,expression_system,protein_sequence,nucleotide_id,nucleotide_sequence,bfactor_or_plddt,secondary_structure\n')
        
    nucleotide_protein_seqs_df.parallel_apply(lambda data: fetch_pdb_and_plddt(data=data, output_path=args.output_path, name=args.name, download=args.download), axis=1)
    
    if not args.download:
        shutil.rmtree(f'{args.output_path}/cif_files/')
        shutil.rmtree(f'{args.output_path}/pdb_files/')
        
if __name__ == "__main__":
    main()