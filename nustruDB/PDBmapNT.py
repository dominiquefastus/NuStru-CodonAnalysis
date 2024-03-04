#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import requests
import logging
import json
import os

import pandas as pd
import multiprocessing as mp
from itertools import islice

from tqdm import tqdm

import mysql.connector
from mysql.connector import Error
from getpass import getpass

from Bio import Entrez
from biopandas.pdb import PandasPdb 

def connect_DB():
    nustruDB = mysql.connector.connect(
        host="localhost",
        user=input("Enter username: "),
        password=getpass("Enter password: "),
        database="nustruDB"
    )
    
    return nustruDB
    
def execute_database(DB, method, table, source, entry_id, gene_name, organism, expression_system, mitochondrial, protein_sequence, nucleotide_id, nucleotide_sequence, plddt):
    if DB is None:
       logging.error("Database connection is not established.")
       exit(1)

    cursor = DB.cursor()

    if method == "INSERT":
        try:
            insert_entry = '''INSERT IGNORE INTO {} 
                            (source, primary_id, gene_name, organism, expression_system, mitochondrial, protein_sequence, nucleotide_id, nucleotide_sequence, plddt) 
                            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'''.format(table)
            entry = (source, entry_id, gene_name, organism, expression_system, mitochondrial, protein_sequence, nucleotide_id, nucleotide_sequence, plddt)
            
            cursor.execute(insert_entry, entry)
            DB.commit()
        except:
            logging.info(f'Entry {entry_id} not successfully inserted in {table} of SQL database!')
            pass
        
    elif method == "UPDATE":
        try:
            update_entry = '''UPDATE {} 
                            SET source = %s, primary_id = %s, gene_name = %s, organism = %s, expression_system = %s,
                            mitochondrial = %s, protein_sequence = %s, nucleotide_id = %s, nucleotide_sequence = %s, 
                            plddt = %s, WHERE primary_id = %s'''.format(table)
            entry = (source, entry_id, gene_name, organism, expression_system, mitochondrial, protein_sequence, nucleotide_id, nucleotide_sequence, plddt)
            
            cursor.execute(update_entry, entry)
            DB.commit()
        except:
           logging.info(f'Entry {entry_id} not successfully updated in {table}!')
           pass
        
def insert_pandas(df, source, entry_id, gene_name, organism, expression_system, mitochondrial, protein_sequence, nucleotide_id, nucleotide_sequence, plddt):
    try:
        df = df._append({"source": source, "primary_id": entry_id, "gene_name": gene_name, "organism": organism, "expression_system": expression_system, "mitochondrial": mitochondrial,
                        "protein_sequence": protein_sequence, "nucleotide_id": nucleotide_id, "nucleotide_sequence": nucleotide_sequence, "plddt": plddt}, ignore_index=True)
    except:
        logging.error(f'Entry {entry_id} not successfully inserted in pandas DataFrame!')
        pass
    
    return df

def get_base_data(entryID):
    try:
        url = rf"""https://data.rcsb.org/graphql?query=%7B%0A%20%20entry(entry_id%3A%20%22{entryID}%22)%20%7B%0A%20%20%20%20rcsb_id%0A%20
        %20%20%20polymer_entities%20%7B%0A%20%20%20%20%20%20rcsb_entity_source_organism%20%7B%0A%20%20%20%20%20%20%20%20scientific_name
        %0A%20%20%20%20%20%20%20%20ncbi_scientific_name%0A%20%20%20%20%20%20%20%20rcsb_gene_name%20%7B%0A%20%20%20%20%20%20%20%20%20
        %20value%0A%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20rcsb_entity_host_organism%20%7B%0A%20%20%20%20
        %20%20%20%20ncbi_scientific_name%0A%20%20%20%20%20%20%7D%0A%20%20%20%20%7D%0A%20%20%7D%0A%7D%0A"""
        
        response = requests.get(url=url)

        if response.status_code == 200:
            response = json.loads(response.text)

            try:
                organism = response['data']['entry']['polymer_entities'][0]['rcsb_entity_source_organism'][0]['ncbi_scientific_name']
            except:
                organism = "NaN"
        
            try:
                gene_name = response['data']['entry']['polymer_entities'][0]['rcsb_entity_source_organism'][0]['rcsb_gene_name'][0]['value']
            except:
                gene_name = "NaN"
                
            try:
                expression_system = response['data']['entry']['polymer_entities'][0]['rcsb_entity_host_organism'][0]['ncbi_scientific_name']
            except:
                expression_system = "NaN"
    except:
        logging.error(f"Base data for protein {entryID} not available!")
        pass
        
    return organism, gene_name, expression_system
       
def get_allignment(entryID, id_type):
    allignment_range = {}
    exon_shift_range = {}
    
    if id_type == "pdb":
        pdbID = entryID
        
        query = '{alignment(from: PDB_ENTITY, to: NCBI_GENOME, queryId: "%s" ) { query_sequence target_alignment { target_id orientation aligned_regions { query_begin query_end target_begin target_end exon_shift } } }}' % pdbID
        
    elif id_type == "uniprot":
        uniprotID = entryID
        
        query = '{alignment(from: UNIPROT, to: NCBI_GENOME, queryId: "%s" ) { query_sequence target_alignment { target_id orientation aligned_regions { query_begin query_end target_begin target_end exon_shift } } }}' % uniprotID
        
    else:
        raise Exception()
        print("Specify provided id")
        exit(1)
        
    # define graphql url from pdb
    url = f'https://1d-coordinates.rcsb.org/graphql?query={query}'
    
    response = requests.get(url=url)
    
    # print("response status code: ", response.status_code)
    
    if response.status_code == 200:
        response = json.loads(response.text)
        
        pdb_sequence = response['data']['alignment']['query_sequence']
        genomeID = response['data']['alignment']['target_alignment'][0]['target_id']
        oritentation = response['data']['alignment']['target_alignment'][0]['orientation']
        
        # loop through the sequence positions
        for id, position in enumerate(response['data']['alignment']['target_alignment'][0]['aligned_regions']):
            # print(position)
            range = [position['target_begin'], position['target_end']]
            exon_shift_range[id] = position['exon_shift']
            allignment_range[id] = range
      
    return pdb_sequence, genomeID, oritentation, allignment_range, exon_shift_range
    
        
def map_uniprot(pdbID):
    query = '{alignment(from: PDB_ENTITY, to: UNIPROT, queryId: "%s" ) { query_sequence target_alignment { target_id target_sequence } }}' % pdbID
    # define graphql url from pdb
    url = f'https://1d-coordinates.rcsb.org/graphql?query={query}'
    
    response = requests.get(url=url)
    
    if response.status_code == 200:
        response = json.loads(response.text)
        uniprot_id = response['data']['alignment']['target_alignment'][0]['target_id']
        
    return uniprot_id
    
def retrieve_nucleotide_seqs(genomeID=None, orientation=None, seqSTART=None, seqEND=None):
    # same as requesting by url
    # f"https://www.ncbi.nlm.nih.gov/nuccore/{genomeID}?report=fasta&from={seqSTART}&to={seqEND}"
    Entrez.email = "dominique.fastus@biochemistry.lu.se"    
    
    if seqSTART > seqEND:
        seqSTART, seqEND = seqEND, seqSTART
        
    # set strand orientation (positive or negative strand)
    if orientation == -1:
        orientation = 2
    else:
        orientation = 1
           
    handle = Entrez.efetch(db="nuccore", id=genomeID, seq_start=seqSTART, seq_stop=seqEND, strand=orientation, rettype="fasta", retmode="text")
    
    return handle.read().partition("\n")[2].strip()

def get_plddt(uniprotID):
    try:
        ppdb = PandasPdb().fetch_pdb(uniprot_id=uniprotID, source="alphafold2-v2")
        residue_b_factors = ppdb.df["ATOM"][["b_factor","residue_number"]]


        residue_b_factors = residue_b_factors.groupby(["residue_number"]).mean().round(4)

        residue_b_factors_dict = residue_b_factors.to_dict()
    except:
        logging.error(f'No plddt value for protein {uniprotID} available')
        residue_b_factors_dict = {"b_factor": "NaN"}

    return residue_b_factors_dict["b_factor"]

def main():
    parser = argparse.ArgumentParser(
        prog="PDBmapNT",
        description="Maps PDB ID to nucleotide sequence and prints an allignment of the pdb protein sequence to the nucleotide sequence"
    )
    parser.add_argument(
        '-i', '--input', required=True, type=str, dest="entryID",
        help="Path to file containing PDB IDs with comma separation."
        )
    parser.add_argument(
        "--sql", action="store_true", dest="sql",
        help="Store the data in a SQL database."
        )
    parser.add_argument(
        "--pandas", dest="pandas", action="store_true",
        help="Store the data in a pandas DataFrame.")
    parser.add_argument( 
        '-o', '--output', type=str, dest="output_path", required=True,
        help='Output path or directory to store the log file and the data.'
    )
    parser.add_argument(
        '-n', '--name', type=str, dest="name", required=True,
        help='Name of the output files and log file.'
    )
    parser.add_argument(
        '--map-uniprot', action="store_true", dest="map_uniprot", default=True,
        help="Map uniprot ID to nucleotide sequence."
    )
    args = parser.parse_args()
    
    if args.sql:
        nustruDB = connect_DB()
    elif args.pandas:
        nucleotide_protein_seqs_df = pd.DataFrame(columns=["source", "primary_id", "gene_name", "organism", "expression_system", "mitochondrial", "protein_sequence", "nucleotide_id", "nucleotide_sequence"])
        if not os.path.exists(args.pandas):
            nucleotide_protein_seqs_df.to_csv(f'{args.output_path}/{args.name}.csv', mode='w', index=False, header=True)
    else:
        raise Exception("Please provide a way to store the data.")
        exit(1)
        
    logging.basicConfig(filename=f'{args.output_path}/{args.name}.log',
                filemode='a',
                format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                datefmt='%H:%M:%S',
                level=logging.ERROR)
        
    with open(args.entryID,'r') as entryIDs_file:
        entryIDs_file = entryIDs_file.read()
        
        progress = tqdm(total=len(entryIDs_file.replace(" ", "").split(',')))

        for pdb_id in entryIDs_file.replace(" ", "").split(','):
            organism, gene_name, expression_system = get_base_data(entryID=pdb_id)
            
            pdb_entry = pdb_id
            pdb_id = pdb_id + "_1"
            
            if args.map_uniprot:
                try:
                    uniprotID = map_uniprot(pdb_id)
                except:
                    logging.error(f'No uniprot ID for protein {pdb_id} available')
                    uniprotID = None
            
            try:
                pdb_sequence, genomeID, oritentation, allignment_range, exon_shift_range = get_allignment(entryID=pdb_id, id_type="pdb")
                
                nu_sequence = ""
                for(seqSTART, seqEND), (exon_range) in zip(allignment_range.values(), exon_shift_range.values()):
                    nu_sequence += retrieve_nucleotide_seqs(genomeID=genomeID, orientation=oritentation, seqSTART=seqSTART, seqEND=seqEND)
                    
                    if len(exon_range) == 2:
                        exonSTART = exon_range[0]
                        exonEND = exon_range[1]
                        
                        nu_sequence += retrieve_nucleotide_seqs(genomeID=genomeID, orientation=oritentation, seqSTART=exonSTART, seqEND=exonEND)
                        
                    elif len(exon_range) == 1:
                        exonSTART = exon_range[0]
                        exonEND = exon_range[0]
                    
                        
                        nu_sequence += retrieve_nucleotide_seqs(genomeID=genomeID, orientation=oritentation, seqSTART=exonSTART, seqEND=exonEND)
                         
                    elif len(exon_range) == 0:
                        pass
                    
                    else:
                        print(exon_range)
                    
                nu_sequence = nu_sequence.replace('\n','')
                
                if args.sql:
                    execute_database(DB=nustruDB, method="INSERT", table="nucleotide_protein_seqs", source="pdb", entry_id=pdb_entry, gene_name=gene_name, organism=organism, 
                                     expression_system=expression_system, mitochondrial="False", protein_sequence=pdb_sequence,
                                     nucleotide_id=genomeID, nucleotide_sequence=nu_sequence, plddt="NaN")
                else:
                    nucleotide_protein_seqs_df = insert_pandas(df=nucleotide_protein_seqs_df, source="pdb", entry_id=pdb_entry, gene_name=gene_name, organism=organism, 
                                                               expression_system=expression_system, mitochondrial="False", protein_sequence=pdb_sequence,
                                                               nucleotide_id=genomeID, nucleotide_sequence=nu_sequence, plddt="NaN")
                    
                    nucleotide_protein_seqs_df.to_csv(f'{args.output_path}/{args.name}.csv', mode='a', index=False, header=False)
                    del nucleotide_protein_seqs_df
                    
                
            except:
                logging.error(f"Genomic coordinates for protein {pdb_id} not available!")
                continue
            
            nucleotide_protein_seqs_df = pd.DataFrame(columns=["source", "primary_id", "gene_name", "organism", "expression_system", "mitochondrial", "protein_sequence", "nucleotide_id", "nucleotide_sequence"])
                
            if args.map_uniprot:
                try:
                    uniprot_sequence, genomeID, oritentation, allignment_range, exon_shift_range = get_allignment(entryID=uniprotID, id_type="uniprot")
                    
                    nu_sequence = ""
                    for(seqSTART, seqEND), (exon_range) in zip(allignment_range.values(), exon_shift_range.values()):
                        nu_sequence += retrieve_nucleotide_seqs(genomeID=genomeID, orientation=oritentation, seqSTART=seqSTART, seqEND=seqEND)
                        
                        if len(exon_range) == 2:
                            exonSTART = exon_range[0]
                            exonEND = exon_range[1]
                            
                            nu_sequence += retrieve_nucleotide_seqs(genomeID=genomeID, orientation=oritentation, seqSTART=exonSTART, seqEND=exonEND)
                            
                        elif len(exon_range) == 1:
                            exonSTART = exon_range[0]
                            exonEND = exon_range[0]
                        
                            
                            nu_sequence += retrieve_nucleotide_seqs(genomeID=genomeID, orientation=oritentation, seqSTART=exonSTART, seqEND=exonEND)
                            
                        elif len(exon_range) == 0:
                            pass
                        
                        else:
                            logging.error(f'The exon range was wrongly parsed, as it is: {exon_range}')
                        
                    nu_sequence = nu_sequence.replace('\n','')
                    
                    plddt = get_plddt(uniprotID=uniprotID)
                    
                    if args.sql:
                        execute_database(DB=nustruDB, method="INSERT", table="nucleotide_protein_seqs", source="uniprot", entry_id=uniprotID, gene_name=gene_name, organism=organism, 
                                        expression_system=expression_system, mitochondrial="False", protein_sequence=uniprot_sequence,
                                        nucleotide_id=genomeID, nucleotide_sequence=nu_sequence, plddt="NaN")
                    else:
                        nucleotide_protein_seqs_df = insert_pandas(df=nucleotide_protein_seqs_df, source="uniprot", entry_id=uniprotID, gene_name=gene_name, organism=organism, 
                                                                expression_system=expression_system, mitochondrial="False", protein_sequence=uniprot_sequence,
                                                                nucleotide_id=genomeID, nucleotide_sequence=nu_sequence, plddt=plddt)
                        
                        nucleotide_protein_seqs_df.to_csv(f'{args.output_path}/{args.name}.csv', mode='a', index=False, header=False)
                        del nucleotide_protein_seqs_df
                        
                except:
                    logging.error(f"Genomic coordinates for protein {uniprotID} not available!")
                    continue
            
            nucleotide_protein_seqs_df = pd.DataFrame(columns=["source", "primary_id", "gene_name", "organism", "expression_system", "mitochondrial", "protein_sequence", "nucleotide_id", "nucleotide_sequence"])
            progress.update(1)            

if __name__ == '__main__':
    main()
