import requests
import json

import asyncio

from multiprocessing import Pool

import mysql.connector
from mysql.connector import Error

import sqlite3
from sqlite3 import Error

from Bio import Entrez

nustruDB = mysql.connector.connect(
    host="localhost",
    user="root",
    database="nustruDB"
)

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
    
def retrieve_nucleotide_seq(ncbiDB="nuccore", genomeID=None, rettype="fasta", retmode="text"):
    Entrez.email = "dominique.fastus@biochemistry.lu.se"    
    handle = Entrez.efetch(db=ncbiDB, id=genomeID, rettype=rettype, retmode=retmode)
    
    return handle.read()


#------------------------------------------------

base_data = 'https://rest.uniprot.org/uniprotkb/search?query=Q9UJV3&fields=gene_primary,sequence,organism_name,cc_subcellular_location'
pdbs_ref = 'https://rest.uniprot.org/uniprotkb/search?query=Q9UJV3&fields=xref_pdb'
url_genomic_locations = 'https://rest.uniprot.org/uniprotkb/search?query=Q9UJV3&fields=xref_embl' 


base_data_response = requests.get(base_data)
if base_data_response.status_code == 200:
    print(base_data_response.text)


response = requests.get(url_genomic_locations)
ref_available = False
    
if response.status_code == 200:
    print(response.text)
    
    ids = []
    response = json.loads(response.text)
    if response['results'][0]['uniProtKBCrossReferences']:
        ref_available = True
        for ccds_ref in response['results'][0]['uniProtKBCrossReferences']:
            ccds_id = ccds_ref['id']
            ids.append(ccds_id)
            try:
                isoformId = ccds_ref['isoformId']
                print(ccds_id, isoformId)
            except:
                print(ccds_id)
    
    
    handle = retrieve_nucleotide_seq(ncbiDB="nuccore",id=ids[0], rettype="fasta_cds_na")
    print(handle)