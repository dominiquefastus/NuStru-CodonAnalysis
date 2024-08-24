#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
INFO: Pandarallel will run on 10 workers.
INFO: Pandarallel will use standard multiprocessing data transfer (pipe) to transfer data between the main process and workers.
usage: db_fetch.py [-h] -i INPUT_FILE -o OUTPUT_PATH -n NAME [-d] [-w]

Fetch the pdb and map the secondary structure, as well as plddt or bfactor to the protein sequence.

options:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input INPUT_FILE
                        Input file of csv formatted uniprot or pdb entries.
  -o OUTPUT_PATH, --output OUTPUT_PATH
                        Output to store the new csv with secondary structure information.
  -n NAME, --name NAME  Name of the output files and log file.
  -d, --download        Download the pdb files.
  -w, --overwrite       If file name already exists, overwrite it. Default is False.
  
Script to retrieve secondary structure information from pdb files and map it to the protein sequence. It also maps the plddt or bfactor to the protein sequence.
Both data are stored in a new csv file with the respective columns.
Example: python db_fetch.py -i Example/examples_nustruDB/example_dbfetch.csv -o . -n example_dbfetch_nustru_fetched
"""

import re
import shutil
import argparse

import logging
import pandas as pd

from pathlib import Path

import requests
from requests.adapters import HTTPAdapter, Retry

from biopandas.pdb import PandasPdb 

from Bio.PDB import MMCIFParser
cifp = MMCIFParser(QUIET=True) # initialize the MMCIFParser and silence the warnings 
    
from Bio.PDB.DSSP import DSSP
from Bio.SeqUtils import seq1

from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True, nb_workers=10, verbose=0) 

class sequenceFeautureFetcher:
    def __init__(self):
        # create a session with retries, with a total of 6 retries, a backoff factor of 0.2 and status forcelist of 502, 503, 504
        self.session = requests.Session()
        self.retries = Retry(total=6, backoff_factor=0.2, status_forcelist=[ 502, 503, 504 ])
        self.session.mount('https://', HTTPAdapter(max_retries=self.retries)) # mount the session with the retries
          
    def reassign_position(self, sequence, residue_position):
        """Match the residue position between the sequence and pdb file"""
        substring = "".join(residue_position.values()) # join the residues to a string
        res = re.search(substring, sequence) # search for the substring in the sequence
        
        # if the substring is found in the sequence (for alphafold it should always be the same sequence as the pdb file)
        if res:
            # get the position range of the substring in the sequence
            position_range = range(res.start(), res.end())
            new_dict = {}
            number_dict = {}

            # reassign the residue position to the new sequence
            # create a new dictionary with the new residue positions and a number dictionary for the old residue positions
            for position_seq, position_dict in zip(position_range, residue_position):
                new_dict[position_seq + 1] = residue_position[position_dict]
                number_dict[position_dict] = position_seq + 1

            return new_dict, number_dict
        else:
            print("Substring not found in the string.")
            return None

    def fetch_features(self, chain_id, sequence, output_path, model_id):
        """Fetch the bfactor and secondary structure from the pdb file"""
        # read the pdb file with PandasPdb
        ppdb = PandasPdb().read_pdb(f"{output_path}/pdb_files/{model_id}.pdb")
        
        # get the residue position and name from the pdb file
        residue_position = ppdb.df["ATOM"][ppdb.df['ATOM']['chain_id'] == chain_id][["residue_number","residue_name"]]
        
        # group atoms by residue number and get the last residue name
        residue_position = residue_position.groupby(["residue_number"])["residue_name"].last().apply(seq1)
        residue_position_dict = residue_position.to_dict() # convert residue position to dictionary

        # reassign the residue position to the sequence as the pdb file can contain missing residues or different residue positions
        # otherwise the sequence starting from 1 would not match the pdb file
        residue_position_dict, number_dict = self.reassign_position(sequence=sequence, residue_position=residue_position_dict)
        
        # get the bfactor for each residue
        residue_b_factor = ppdb.df["ATOM"][ppdb.df['ATOM']['chain_id'] == chain_id][["b_factor","residue_number"]]
        # group atoms by residue number and get the mean bfactor or plddt (alphafol) for each residue
        residue_b_factor = residue_b_factor.groupby(["residue_number"]).mean().round(4)['b_factor']
        residue_b_factor_dict = residue_b_factor.to_dict() # convert residue bfactor to dictionary
        
        # reassign the residue bfactors to the new residue positions adapted from the pdb file
        # prevents missinterpratetion of the bfactor values (especially for experimental structures)
        updated_residue_b_factor_dict = {}
        for position in residue_b_factor_dict.keys():
            new_position = number_dict[position]
            updated_residue_b_factor_dict[new_position] = residue_b_factor[position]
        
        # get the secondary structure for each residue by parsing the cif file
        # cif file works more reliable than the pdb file for the dssp algorithm
        structure = cifp.get_structure(model_id, f"{output_path}/cif_files/{model_id}.cif")

        # get the secondary structure for each residue
        model = structure[0] # only one model (chain) can be applied at a time
        
        # run the dssp algorithm on the model (has to be installed on the system)
        dssp = DSSP(model, f"{output_path}/cif_files/{model_id}.cif")
        # get the secondary structure for each residue and assign it to the residue positions present in the pdb file
        # normally the pdb contains the same residues as the cif file, but for experimental structures it can differ
        secondary_structure = [dssp_data[0:3:2] for dssp_data in dssp if dssp_data[0] in residue_position_dict.keys()]
        secondary_structure_dict = dict(secondary_structure)
        
        return updated_residue_b_factor_dict, secondary_structure_dict

    def fetch_pdb_and_plddt(self, data, output_path, name, download=False):
        """Fetch the pdb file and plddt or bfactor for the protein sequence"""
        # make temporary directories to store the cif and pdb files
        Path(f'{output_path}/pdb_files').mkdir(parents=True, exist_ok=True)
        Path(f'{output_path}/cif_files').mkdir(parents=True, exist_ok=True)

        try:
            # check if the source is uniprot or pdb
            # for uniprot the model id is the primary id and the chain id is A
            # the sequence should be the same as the alphafold model, so no reassigning to the residue positions needed
            if data['source'] == 'uniprot':
                model_id = data['primary_id']
                chain_id = 'A'
                
                try:
                    # use the alphafold api to retrieve the cif and pdb files
                    # the alphafold versions of models can be different, so the version has to be checked in the future
                    af2_version = 4
                    response_cif = self.session.get(f"https://alphafold.ebi.ac.uk/files/AF-{model_id}-F1-model_v{af2_version}.cif")
                    response_pdb = self.session.get(f"https://alphafold.ebi.ac.uk/files/AF-{model_id}-F1-model_v{af2_version}.pdb")
                    
                    # write the the cif and pdb files to the cif_files and pdb_files directories if the response is not an error
                    if "Error" not in response_cif.text and "Error" not in response_pdb.text:           
                        with open(f"{output_path}/cif_files/{model_id}.cif", "wb") as cif_file, open(f"{output_path}/pdb_files/{model_id}.pdb", "wb") as pdb_file:
                            cif_file.write(response_cif.content)
                            pdb_file.write(response_pdb.content)
        
                        # get the bfactor or plddt and secondary structure for the protein sequence as a dictionary
                        residue_b_factor_dict,secondary_structure_dict = self.fetch_features(model_id=model_id, chain_id=chain_id, sequence=data['protein_sequence'], output_path=output_path)
                        
                    else:
                        logging.error(f"Error: could not retrieve cif or file for {model_id}.")
                        pass
                except: 
                    logging.error(f"Could not assign plddt to {model_id}.")
                    pass
                
            # for pdb each chain of the pdb file needs to be individually processed 
            else:
                # assign the model id and chain id from the primary id
                model_id = data['primary_id'].replace('"','').split("_")[0]
                chain_id = data['primary_id'].replace('"','').split("_")[1]
                chain_id = list(chain_id)[0] # only one chain id is mapped if multiple chains have the same protein sequence
                
                try:
                    # retrieve the cif and pdb files from the rcsb database api
                    response_cif = self.session.get(f"https://files.rcsb.org/download/{model_id}.cif")
                    response_pdb = self.session.get(f"https://files.rcsb.org/download/{model_id}.pdb")
                    
                    # write the the cif and pdb files to the cif_files and pdb_files directories if the response is not an error
                    if "Error" not in response_cif.text and "Error" not in response_pdb.text:        
                        with open(f"{output_path}/cif_files/{model_id}.cif", "wb") as cif_file, open(f"{output_path}/pdb_files/{model_id}.pdb", "wb") as pdb_file:
                            cif_file.write(response_cif.content)
                            pdb_file.write(response_pdb.content)

                        # again get the bfactor or plddt and secondary structure for the protein sequence as a dictionary based on the chain id
                        residue_b_factor_dict,secondary_structure_dict = self.fetch_features(model_id=model_id, chain_id=chain_id, sequence=data['protein_sequence'], output_path=output_path)
                        
                    else:
                        logging.error(f"Error: could not retrieve cif or pdb file for {model_id}.")
                        pass
                except:
                    logging.error(f"Could not assign plddt to {model_id}.")
                    pass
            
            # assign the bfactor or plddt and secondary structure columns of the new dataframe to the data dictionaries
            data['bfactor_or_plddt'] = residue_b_factor_dict
            data['secondary_structure'] = secondary_structure_dict

            # write the new data to the csv file and delete the data and new data to free memory
            new_data = pd.DataFrame([{'source': data['source'], 'primary_id': data['primary_id'], 'gene_name': data['gene_name'], 'organism': data['organism'], 'expression_system': data['expression_system'],
                    'protein_sequence': data['protein_sequence'], 'nucleotide_id': data['nucleotide_id'], 'nucleotide_sequence': data['nucleotide_sequence'], 'bfactor_or_plddt': data['bfactor_or_plddt'], 'secondary_structure': data['secondary_structure']}])
            new_data.to_csv(f'{output_path}/{name}.csv', mode='a', index=False, header=False)
            del data, new_data
            
            # check if the download flag is set to False, if so delete the cif and pdb files
            if not download:
                Path(f'{output_path}/pdb_files/{model_id}.pdb').unlink()
                Path(f'{output_path}/cif_files/{model_id}.cif').unlink()
                
        except:
            logging.error(f"Error: {data['primary_id']} not found.")
            pass

def main():
    parser = argparse.ArgumentParser(
        prog='db_fetch.py',
        description="Fetch the pdb and map the secondary structure, as well as plddt or bfactor to the protein sequence."
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
        help='Download the structure files.'
    )
    parser.add_argument(
        '-w', '--overwrite', action="store_true", dest="overwrite", required=False, default=False,
        help='If file name already exists, overwrite it. Default is False.' 
    )
    parser.add_argument(
        '--delimiter', type=str, dest="delimiter", required=False, default=';',
    )
    args = parser.parse_args()
    
    logging.basicConfig(filename=f'{args.output_path}/{args.name}.log',
            filemode='a',
            format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
            datefmt='%H:%M:%S',
            level=logging.ERROR)
    
    # read the input file
    nucleotide_protein_seqs_df = pd.read_csv(args.input_file, index_col=False, delimiter=args.delimiter)
    print(nucleotide_protein_seqs_df.head())
    
    # create the new csv file with the respective columns if it does not exist or the overwrite flag is set
    if not Path(f'{args.output_path}/{args.name}.csv').exists() or args.overwrite:
        with open(f'{args.output_path}/{args.name}.csv', mode='w') as f:
            f.write('source,primary_id,gene_name,organism,expression_system,protein_sequence,nucleotide_id,nucleotide_sequence,bfactor_or_plddt,secondary_structure\n')
    else:
        print(f"Error: {args.output_path} already exists. Use -w to overwrite.")
        exit(1)
    
    sequenceFF = sequenceFeautureFetcher() # initialize the sequence feature fetcher    
    # apply the fetch_pdb_and_plddt function to get the features to the dataframe in parallel
    nucleotide_protein_seqs_df.parallel_apply(lambda data: sequenceFF.fetch_pdb_and_plddt(data=data, output_path=args.output_path, name=args.name, download=args.download), axis=1)
    
    # delete the cif and pdb directories if the download flag is set to False
    if not args.download:
        shutil.rmtree(f'{args.output_path}/cif_files/')
        shutil.rmtree(f'{args.output_path}/pdb_files/')
        
if __name__ == "__main__":
    main()