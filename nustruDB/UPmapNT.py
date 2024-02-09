#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import requests
import json
import os

from typing import List

import pandas as pd
import multiprocessing as mp
from itertools import islice

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
        print("Error! Database connection is not established.")
        return

    cursor = DB.cursor()
    
    if method == "INSERT":
        insert_entry = '''INSERT IGNORE INTO {} 
                          (source, primary_id, gene_name, organism, expression_system, mitochondrial, protein_sequence, nucleotide_id, nucleotide_sequence, plddt) 
                          VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'''.format(table)
        entry = (source, entry_id, gene_name, organism, expression_system, mitochondrial, protein_sequence, nucleotide_id, nucleotide_sequence, plddt)
        
        cursor.execute(insert_entry, entry)
        DB.commit()
        print(f'Entry {entry_id} successfully inserted in {table} of SQL database!')
        
    elif method == "UPDATE":
        update_entry = '''UPDATE {} 
                          SET source = %s, primary_id = %s, gene_name = %s, organism = %s, expression_system = %s,
                          mitochondrial = %s, protein_sequence = %s, nucleotide_id = %s, nucleotide_sequence = %s, 
                          plddt = %s, WHERE primary_id = %s'''.format(table)
        entry = (source, entry_id, gene_name, organism, expression_system, mitochondrial, protein_sequence, nucleotide_id, nucleotide_sequence, plddt)
        
        cursor.execute(update_entry, entry)
        DB.commit()
        print(f'Entry {entry_id} successfully updated in {table}!')
        
    elif method == "SELECT ALL":
        cursor.execute(f'SELECT * FROM {table}')
        
        rows = cursor.fetchall()

        for row in rows:
            print(row)
    
def insert_pandas(df, source, entry_id, gene_name, organism, expression_system, mitochondrial, protein_sequence, nucleotide_id, nucleotide_sequence, plddt):
    df = df._append({"source": source, "primary_id": entry_id, "gene_name": gene_name, "organism": organism, "expression_system": expression_system, "mitochondrial": mitochondrial, 
                    "protein_sequence": protein_sequence, "nucleotide_id": nucleotide_id, "nucleotide_sequence": nucleotide_sequence, "plddt": plddt}, ignore_index=True)
    
    print(f'Entry {entry_id} successfully inserted in pandas DataFrame!')
    
    return df

def retrieve_nucleotide_seq(ncbiDB="nuccore", entryID=None, rettype="fasta", retmode="text"):
    Entrez.email = "dominique.fastus@biochemistry.lu.se"    
    handle = Entrez.efetch(db=ncbiDB, id=entryID, rettype=rettype, retmode=retmode)
    
    return handle.read()

pdbs_ref = 'https://rest.uniprot.org/uniprotkb/search?query=Q9UJV3&fields=xref_pdb'

def get_base_data(uniprotID):
    base_data = 'https://rest.uniprot.org/uniprotkb/search?query=%s&fields=gene_primary,sequence,organism_name,cc_subcellular_location' %uniprotID
    
    response = requests.get(base_data)
    # print(response)
        
    if response.status_code == 200:
        
        response = json.loads(response.text)
        
        organism = response['results'][0]['organism']['scientificName']
        gene_name = response['results'][0]['genes'][0]['geneName']['value']
        
        sequence = response['results'][0]['sequence']['value']
    
    return organism,gene_name,sequence

def filter_sequence(sequences: str, searchID: List) -> None:
    search_criteria = searchID
    unique_sequences = set()
    matching_data = []
    entries = sequences.split('>')

    for entry in entries:
        if any(criterion in entry for criterion in search_criteria):
            header, sequence = entry.split('\n', 1)
            sequence = sequence.strip()
            # Check if the sequence is already added to ensure uniqueness
            if sequence not in unique_sequences:
                unique_sequences.add(sequence)
                matching_data.append((header, sequence.replace('\n',"")))

    return matching_data

def check_isoform(uniprotID):
    url_gene = 'https://rest.uniprot.org/uniprotkb/search?query=%s&fields=cc_alternative_products' %uniprotID

    response_gene = requests.get(url_gene)

    if response_gene.status_code == 200: 
        response = json.loads(response_gene.text)
        
        try:
            if any(response['results'][0]['comments'][0]['isoforms']):
                return True
        except:
            return False

def get_isoform(uniprotID):
    url = 'https://rest.uniprot.org/uniprotkb/search?query=%s&fields=xref_refseq' %uniprotID

    response = requests.get(url=url)

    print("response status code: ", response.status_code)

    if response.status_code == 200:
        response = json.loads(response.text)
        print(response)

        if response['results'][0]['uniProtKBCrossReferences']:
            protein_ids = []
            cds_ids = []
            isoformIds = []
            
            for cds_ref in response['results'][0]['uniProtKBCrossReferences']:
                protein_id = cds_ref['id']
                
                if protein_id.startswith("NP"):
                    protein_ids.append(protein_id)
                    
                    cds_id = cds_ref['properties'][0]['value']

                    if cds_id.startswith("NM"):
                        cds_ids.append(cds_id)
                    
                try:
                    isoformId = cds_ref['isoformId']
                    isoformIds.append(isoformId)
                except:
                    pass
                    
            return protein_ids, cds_ids, isoformIds

def get_cds(uniprotID):
    url_gene = 'https://rest.uniprot.org/uniprotkb/search?query=%s&fields=xref_embl' %uniprotID

    response_gene = requests.get(url_gene)

    if response_gene.status_code == 200: 
        cds_ids = []
        protein_ids = []
        
        response = json.loads(response_gene.text)

        
        if response['results'][0]['uniProtKBCrossReferences']:
            for cds_ref in response['results'][0]['uniProtKBCrossReferences']:
                cds_id = cds_ref['id']
                
                for properties in cds_ref['properties']:

                    if properties['key'] == 'ProteinId':
                        protein_id = properties['value']
                
                    if properties['key'] == 'Status':
                        status_code = properties['value']
                    
                if status_code == "-":
                    cds_ids.append(cds_id)
                    protein_ids.append(protein_id)
                        
        return cds_ids, protein_ids
    
def get_plddt(uniprotID):
    ppdb = PandasPdb().fetch_pdb(uniprot_id=uniprotID, source="alphafold2-v2")
    residue_b_factors = ppdb.df["ATOM"][["b_factor","residue_number"]]


    residue_b_factors = residue_b_factors.groupby(["residue_number"]).mean().round(4)

    residue_b_factors_dict = residue_b_factors.to_dict()

    return residue_b_factors_dict["b_factor"]
    
    
def main():
    parser = argparse.ArgumentParser(
        prog="UPmapNT",
        description="Maps uniprot ID to nucleotide sequence and prints an allignment of the pdb protein sequence to the nucleotide sequence"
    )
    parser.add_argument('entryID')
    parser.add_argument("--sql", action="store_true", dest="sql")
    parser.add_argument("--pandas", dest="pandas")
    args = parser.parse_args()
    
    if args.sql:
        nustruDB = connect_DB()
    elif args.pandas:
        nucleotide_protein_seqs_df = pd.DataFrame(columns=["source", "primary_id", "gene_name", "organism", "expression_system", "mitochondrial", "protein_sequence", "nucleotide_id", "nucleotide_sequence"])
        if os.path.exists(args.pandas):
            nucleotide_protein_seqs_df.to_csv(args.pandas, mode='w', index=False, header=True)
    else:
        print("Please provide a way to store the data.")
        exit(1)
        
    with open(args.entryID,'r') as entryIDs_file:
        entryIDs_file = entryIDs_file.read()

        for uniprotID in entryIDs_file.replace("\n", ",").split(','):
            print(uniprotID)
            
            organism, gene_name, sequence = get_base_data(uniprotID)
    
            try:
                if check_isoform(uniprotID):
                    protein_ids, cds_ids, isoform_ids = get_isoform(uniprotID)
                    
                    isoform_protein_sequences = []
                    for cds_id, protein_id, isoform_id in zip(cds_ids, protein_ids, isoform_ids):
                        print(cds_id, protein_id, isoform_id)
                        nt_response_sequences = retrieve_nucleotide_seq(entryID=cds_id, rettype="fasta_cds_na")
                        
                        nt_response_sequences = filter_sequence(sequences=nt_response_sequences,searchID=[cds_id])
                        
                        protein_response_sequence = retrieve_nucleotide_seq(ncbiDB="protein",entryID=protein_id, rettype="fasta")
                        
                        protein_response_sequence = filter_sequence(sequences=protein_response_sequence,searchID=[protein_id])
                        isoform_protein_sequences.append(protein_response_sequence[0][1])
                        
                        plddt = get_plddt(uniprotID=uniprotID)
                        
                        if args.sql:
                            execute_database(DB=nustruDB, method="INSERT", table="nucleotide_protein_seqs", source="uniprot", entry_id=isoform_id, gene_name=gene_name, organism=organism, 
                                            expression_system="NaN", mitochondrial="False", protein_sequence=sequence,
                                            nucleotide_id=cds_id, nucleotide_sequence=nt_response_sequences[0][1], plddt=plddt)
                        else:
                            nucleotide_protein_seqs_df = insert_pandas(df=nucleotide_protein_seqs_df, source="uniprot", entry_id=isoform_id, gene_name=gene_name, organism=organism, 
                                            expression_system="NaN", mitochondrial="False", protein_sequence=sequence,
                                            nucleotide_id=cds_id, nucleotide_sequence=nt_response_sequences[0][1], plddt=plddt)
                            
                            nucleotide_protein_seqs_df.to_csv(args.pandas, mode='a', index=False, header=False)
                        
                        
                if not check_isoform(uniprotID) or not sequence in isoform_protein_sequences:
                    cds_ids, protein_ids = get_cds(uniprotID)
                    plddt = get_plddt(uniprotID=uniprotID)
                    
                    for cds_id, protein_id in zip(cds_ids, protein_ids):
                        if retrieve_nucleotide_seq(ncbiDB="protein", entryID=protein_id).partition("\n")[2].strip().replace('\n','') == sequence:
                            print(cds_id, protein_id)
                            response_sequences = retrieve_nucleotide_seq(entryID=cds_id, rettype="fasta_cds_na")
                            
                            matched_id, matched_seq = filter_sequence(sequences=response_sequences,searchID=[uniprotID, protein_id])[0]
                            
                            print(matched_seq)
                            print(sequence)
                            
                            if args.sql:
                                execute_database(DB=nustruDB, method="INSERT", table="nucleotide_protein_seqs", source="uniprot", entry_id=uniprotID, gene_name=gene_name, organism=organism, 
                                                expression_system="NaN", mitochondrial="False", protein_sequence=sequence,
                                                nucleotide_id=cds_id, nucleotide_sequence=matched_seq, plddt=plddt)
                            else:
                                nucleotide_protein_seqs_df = insert_pandas(df=nucleotide_protein_seqs_df, source="uniprot", entry_id=uniprotID, gene_name=gene_name, organism=organism, 
                                                expression_system="NaN", mitochondrial="False", protein_sequence=sequence,
                                                nucleotide_id=cds_id, nucleotide_sequence=matched_seq, plddt=plddt)
                                
                                nucleotide_protein_seqs_df.to_csv(args.pandas, mode='a', index=False, header=False)
                                
                            break

            except:
                print(f"Nucleotide sequence for protein {uniprotID} not available!")
                continue
        
if __name__ == '__main__':
    main()