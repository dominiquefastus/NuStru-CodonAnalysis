#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import requests
import json

import mysql.connector
from mysql.connector import Error
from getpass import getpass

from Bio import Entrez

nustruDB = mysql.connector.connect(
    host="localhost",
    user=input("Enter username: "),
    password=getpass("Enter password: "),
    database="nustruDB"
)
    
def execute_database(method, table, source, entry_id, gene_name, organism, expression_system, mitochondrial, protein_sequence, nucleotide_id, nucleotide_sequence):
    if nustruDB is None:
        print("Error! Database connection is not established.")
        return

    cursor = nustruDB.cursor()
    
    if method == "INSERT":
        insert_entry = '''INSERT INTO {} 
                          (source, primary_id, gene_name, organism, expression_system, mitochondrial, protein_sequence, nucleotide_id, nucleotide_sequence) 
                          VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)'''.format(table)
        entry = (source, entry_id, gene_name, organism, expression_system, mitochondrial, protein_sequence, nucleotide_id, nucleotide_sequence)
        
        cursor.execute(insert_entry, entry)
        nustruDB.commit()
        print(f'Entry {entry_id} successfully inserted in {table}!')
        
    elif method == "UPDATE":
        update_entry = '''UPDATE {} 
                          SET protein_sequence = %s, nucleotide_id = %s, nucleotide_sequence = %s 
                          WHERE primary_id = %s'''.format(table)
        entry = (protein_sequence, nucleotide_id, nucleotide_sequence, entry_id)
        
        cursor.execute(update_entry, entry)
        nustruDB.commit()
        print(f'Entry {entry_id} successfully updated in {table}!')
        
    elif method == "SELECT ALL":
        # Select all rows from the table
        cursor.execute(f'SELECT * FROM {table}')
        
        # Fetch all rows from the cursor
        rows = cursor.fetchall()

        # Print the rows
        for row in rows:
            print(row)


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
    
    print("response status code: ", response.status_code)
    
    if response.status_code == 200:
        response = json.loads(response.text)
        
        pdb_sequence = response['data']['alignment']['query_sequence']
        genomeID = response['data']['alignment']['target_alignment'][0]['target_id']
        oritentation = response['data']['alignment']['target_alignment'][0]['orientation']
        
        # loop through the sequence positions
        for id, position in enumerate(response['data']['alignment']['target_alignment'][0]['aligned_regions']):
            print(position)
            range = [position['target_begin'], position['target_end']]
            exon_shift_range[id] = position['exon_shift']
            allignment_range[id] = range
      
    return pdb_sequence, genomeID, oritentation, allignment_range, exon_shift_range
    
        
def map_uniprot(pdbID):
    query = '{alignment(from: PDB_ENTITY, to: UNIPROT, queryId: "%s" ) { query_sequence target_alignment { target_id target_sequence } }}' % pdbID
    # define graphql url from pdb
    url = f'https://1d-coordinates.rcsb.org/graphql?query={query}'
    
    response = requests.get(url=url)
    
    print("response status code: ", response.status_code)
    
    if response.status_code == 200:
        response = json.loads(response.text)
        uniprot_id = response['data']['alignment']['target_alignment'][0]['target_id']
        
    return uniprot_id
    
def retrieve_nucleotide_seqs(genomeID=None, seqSTART=None, seqEND=None, orientation=None):
    # same as requesting by url
    # f"https://www.ncbi.nlm.nih.gov/nuccore/{genomeID}?report=fasta&from={seqSTART}&to={seqEND}"
    
    # set strand orientation (positive or negative strand)
    if orientation == -1:
        orientation = 2
    else:
        orientation = 1
        
    Entrez.email = "dominique.fastus@biochemistry.lu.se"    
    handle = Entrez.efetch(db="nuccore", id=genomeID, seq_start=seqSTART, seq_stop=seqEND, strand=orientation, rettype="fasta", retmode="text")
    
    return handle.read().partition("\n")[2].strip()

def main():
    parser = argparse.ArgumentParser(
        prog="PDBmapNT",
        description="Maps PDB ID to nucleotide sequence and prints an allignment of the pdb protein sequence to the nucleotide sequence"
    )
    parser.add_argument('entryID')
    args = parser.parse_args()
    
    with open(args.entryID,'r') as entryIDs_file:
        entryIDs_file = entryIDs_file.read()

        for pdb_id in entryIDs_file.split(','):
            pdb_id = pdb_id + "_1"
            
            try:
                uniprotID = map_uniprot(pdb_id)
            except:
                print(f'No uniprot ID for protein {pdb_id} available')
                uniprotID = None
            
            try:
                pdb_sequence, genomeID, oritentation, allignment_range, exon_shift_range = get_allignment(entryID=pdb_id, id_type="pdb")
                
                print(pdb_sequence, genomeID, oritentation, allignment_range, exon_shift_range)
                
                nu_sequence = ""
                for(seqSTART, seqEND), (exon_range) in zip(allignment_range.values(), exon_shift_range.values()):
                    nu_sequence += retrieve_nucleotide_seqs(genomeID=genomeID, seqSTART=seqSTART, seqEND=seqEND)
                    
                    if len(exon_range) == 2:
                        exonSTART = exon_range[0]
                        exonEND = exon_range[1]
                        
                        nu_sequence += retrieve_nucleotide_seqs(genomeID=genomeID, seqSTART=exonSTART, seqEND=exonEND)
                        
                    elif len(exon_range) == 1:
                        exonSTART = exon_range[0]
                        exonEND = exon_range[0]
                    
                        
                        nu_sequence += retrieve_nucleotide_seqs(genomeID=genomeID, seqSTART=exonSTART, seqEND=exonEND)
                        
                    elif len(exon_range) == 0:
                        pass
                    
                    else:
                        print(exon_range)
                    
                nu_sequence = nu_sequence.replace('\n','')
                print(nu_sequence, '\n')
                
                execute_database(method="INSERT", table="nucleotide_protein_seqs", source="pdb", entry_id=pdb_id, gene_name="NaN", organism="NaN", 
                                 expression_system="NaN", mitochondrial="False", protein_sequence=pdb_sequence,
                                 nucleotide_id=genomeID, nucleotide_sequence=nu_sequence)
            except Error as e:
                print(e)
                print(f"Genomic coordinated for protein {pdb_id} not available!")
            
            
            print("\n\n")
                
            
            try:
                pdb_sequence, genomeID, oritentation, allignment_range, exon_shift_range = get_allignment(entryID=uniprotID, id_type="uniprot")
                
                print(pdb_sequence, genomeID, oritentation, allignment_range, exon_shift_range)
                
                nu_sequence = ""
                for(seqSTART, seqEND), (exon_range) in zip(allignment_range.values(), exon_shift_range.values()):
                    nu_sequence += retrieve_nucleotide_seqs(genomeID=genomeID, seqSTART=seqSTART, seqEND=seqEND)
                    
                    if len(exon_range) == 2:
                        exonSTART = exon_range[0]
                        exonEND = exon_range[1]
                        
                        nu_sequence += retrieve_nucleotide_seqs(genomeID=genomeID, seqSTART=exonSTART, seqEND=exonEND)
                        
                    elif len(exon_range) == 1:
                        exonSTART = exon_range[0]
                        exonEND = exon_range[0]
                    
                        
                        nu_sequence += retrieve_nucleotide_seqs(genomeID=genomeID, seqSTART=exonSTART, seqEND=exonEND)
                        
                    elif len(exon_range) == 0:
                        pass
                    
                    else:
                        print(exon_range)
                    
                nu_sequence = nu_sequence.replace('\n','')
                print(nu_sequence, '\n')
                
                execute_database(method="INSERT", table="nucleotide_protein_seqs", source="uniprot", entry_id=uniprotID, gene_name="NaN", organism="NaN", 
                                 expression_system="NaN", mitochondrial="False", protein_sequence=pdb_sequence,
                                 nucleotide_id=genomeID, nucleotide_sequence=nu_sequence)
            except Error as e:
                print(e)
                print(f"Genomic coordinated for protein {uniprotID} not available!")
                pass
    

if __name__ == '__main__':
    main()