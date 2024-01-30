import argparse
import requests
import json

import mysql.connector
from mysql.connector import Error

import sqlite3
from sqlite3 import Error

from Bio import Entrez
from Bio import SeqIO

"""
nustruDB = mysql.connector.connect(
    host="localhost",
    user="root",
    passwd="AL-16-C988",
    database="nustruDB"
)
"""

def create_connection(db_file):
    """
    Create database connection to nustruDB
    """
    
    sqliteConnection = None
    try:
        sqliteConnection = sqlite3.connect(db_file)
        sqliteCursor = sqliteConnection.cursor()

        # Create a table
        sqliteCursor.execute('''CREATE TABLE IF NOT EXISTS pdb_entries
                        (primary_id TEXT PRIMARY KEY,
                        secondary_id TEXT SECONDARY KEY, 
                        gene_name TEXT,
                        organism TEXT,
                        expression_system TEXT,
                        mitochondrial TEXT,
                        protein_sequence TEXT,
                        genome_id TEXT, 
                        nucleotide_sequence TEXT);''')
        
        # Close the connection
        sqliteConnection.close()
    
    except Error as e:
        print('Error occurred - ', e)
    finally:
        if sqliteConnection:
            sqliteConnection.close()
            print('SQLite Connection closed')
            
def execute_database(method, db_file, table, entry_id, sec_id, protein_sequence, genome_id, nucleotide_sequence):
    sqliteConnection = sqlite3.connect(db_file)
    sqliteCursor = sqliteConnection.cursor()
    
    if method == "INSERT":
        sqliteCursor.execute(f'INSERT OR IGNORE INTO {table} VALUES (?,?,?,?,?)', (entry_id, sec_id, protein_sequence, genome_id, nucleotide_sequence))
        print(f'Entry {entry_id} succesfull inserted to {table} in nustru.db!')
        
    elif method == "UPDATE":
        sqliteCursor.execute(f'UPDATE {table} SET protein_sequence = ?, genome_id = ?, nucleotide_sequence = ? WHERE pdb_id = ?', (protein_sequence, genome_id, nucleotide_sequence, entry_id))
        
    elif method == "SELECT ALL":
        # Select all rows from the table
        sqliteCursor.execute(f'SELECT * FROM {table}')
        
        # Fetch all rows from the cursor
        rows = sqliteCursor.fetchall()

        # Print the rows
        for row in rows:
            print(row)
    
    sqliteConnection.commit() 
    sqliteConnection.close()
    
def retrieve_nucleotide_seq(ncbiDB="nuccore", entryID=None, rettype="fasta", retmode="text"):
    Entrez.email = "dominique.fastus@biochemistry.lu.se"    
    handle = Entrez.efetch(db=ncbiDB, id=entryID, rettype=rettype, retmode=retmode)
    
    return handle.read()


#------------------------------------------------


pdbs_ref = 'https://rest.uniprot.org/uniprotkb/search?query=Q9UJV3&fields=xref_pdb'

def get_base_data(uniprotID):
    base_data = 'https://rest.uniprot.org/uniprotkb/search?query=%s&fields=gene_primary,sequence,organism_name,cc_subcellular_location' %uniprotID
    
    response = requests.get(base_data)
        
    if response.status_code == 200:
        print(response.text)
        
        response = json.loads(response.text)
        
        organism = response['results'][0]['organism']['scientificName']
        gene_name = response['results'][0]['genes'][0]['geneName']['value']
        
        sequence = response['results'][0]['sequence']['value']
        print(organism,gene_name,sequence)
    
def get_cds(uniprotID):
    url_ccds = 'https://rest.uniprot.org/uniprotkb/search?query=%s&fields=xref_ccds' %uniprotID
    url_gene = 'https://rest.uniprot.org/uniprotkb/search?query=%s&fields=xref_embl' %uniprotID

    response_ccds = requests.get(url_ccds)
      
    if response_ccds.status_code == 200:
        response = json.loads(response_ccds.text)
        
        if response['results'][0]['uniProtKBCrossReferences']:
            for ccds_ref in response['results'][0]['uniProtKBCrossReferences']:
                ccds_id = ccds_ref['id']
                try:
                    isoformId = ccds_ref['isoformId']
                    print(ccds_id, isoformId)
                except:
                    print(ccds_id)
                
        else:
            response_gene = requests.get(url_gene)
        
            if response_gene.status_code == 200: 
                cds_ids = []
                protein_ids = []
                
                response = json.loads(response_gene.text)
                print(response)
                
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
    
    
def main():
    parser = argparse.ArgumentParser(
        prog="PDBmapNT",
        description="Maps Uniprot ID to nucleotide sequence and prints an allignment of the pdb protein sequence to the nucleotide sequence"
    )
    parser.add_argument('entryID')
    args = parser.parse_args()
    
    create_connection("nustru.db")
    
    with open(args.entryID,'r') as entryIDs_file:
        entryIDs_file = entryIDs_file.read()

        for uniprot_id in entryIDs_file.split(','):
            
            
            handle = retrieve_nucleotide_seq(ncbiDB="nuccore",entryID=ids[4], rettype="fasta_cds_na")
            print(handle)
    
if __name__ == '__main__':
    get_base_data("P10486")
    cds_ids, protein_ids = get_cds("P02185")

    
    print(retrieve_nucleotide_seq(entryID=cds_ids, rettype="fasta_cds_na"))