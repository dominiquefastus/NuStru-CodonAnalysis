#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script is slow and not optimized for large datasets. However it provides isoform data mapping.
For fasta uniprot data retrieval, use fastUPmapNT.py instead.

usage: UPmapNT [-h] -i ENTRYID [--sql] [--pandas] -o OUTPUT_PATH -n NAME [-w] -m API_MAIL -k API_KEY

Map uniprot ID to nucleotide sequence and prints an allignment of the pdb protein sequence to the nucleotide sequence

options:
  -h, --help            show this help message and exit
  -i ENTRYID, --input ENTRYID
                        Path to file containing PDB IDs with comma separation.
  --sql                 Store the data in a SQL database.
  --pandas              Store the data in a pandas DataFrame.
  -o OUTPUT_PATH, --output OUTPUT_PATH
                        Output path or directory to store the log file and the data.
  -n NAME, --name NAME  Name of the output files and log file.
  -w, --overwrite       If file name already exists, overwrite it. Default is False.
  -m API_MAIL, --mail API_MAIL
                        Provide mail from ncbi api account: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
  -k API_KEY, --key API_KEY
                        Provide api key from ncbi api account: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/

The script will create a csv file where the data from the uniprot ID mapping is stored. 
The script will also create fasta files with the nucleotide and protein sequences of the PDB and uniprot if the respective flags are set.
Example: python UPmapNT.py -i Example/examples_nustruDB/example_uniprotIDs.txt --pandas -o . -n example_uniprotIDs_nustru --create-fasta all -w
"""

import argparse
import requests
import logging
import json

from typing import List

import pandas as pd

from tqdm import tqdm
from pathlib import Path

import mysql.connector
from getpass import getpass
from requests.adapters import HTTPAdapter, Retry

# create a session with retries, with a total of 6 retries, a backoff factor of 0.2, and status codes 502, 503, and 504
session = requests.Session()
retries = Retry(total=6, backoff_factor=0.2, status_forcelist=[ 502, 503, 504 ])
session.mount('https://', HTTPAdapter(max_retries=retries))

from Bio import Entrez

def connect_DB():
    """Connects to the SQL database."""
    # /* might be removed or extended in future */
    # set up the connection to the database which is hosted locally
    nustruDB = mysql.connector.connect(
        host="localhost",
        user=input("Enter username: "),
        password=getpass("Enter password: "),
        database="nustruDB"
    )
    
    return nustruDB

    
def execute_database(DB, method, table, source, entry_id, gene_name, organism, expression_system, mitochondrial, protein_sequence, nucleotide_id, nucleotide_sequence):
    """Inserts or updates the data in the SQL database."""
    # /* might be removed or extended in future */
    if DB is None:
        logging.error("Error! Database connection is not established.")
        return

    cursor = DB.cursor() # create a cursor object to execute SQL queries
    
    if method == "INSERT":
        try:
            # insert the data into the table and ignore the entry if it already exists
            insert_entry = '''INSERT IGNORE INTO {} 
                            (source, primary_id, gene_name, organism, expression_system, mitochondrial, protein_sequence, nucleotide_id, nucleotide_sequence) 
                            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'''.format(table)
            entry = (source, entry_id, gene_name, organism, expression_system, mitochondrial, protein_sequence, nucleotide_id, nucleotide_sequence)
            
            # execute the query and commit the changes to the database
            cursor.execute(insert_entry, entry)
            DB.commit()
        except:
            logging.error(f'Entry {entry_id} not successfully inserted in {table} of SQL database!')
            pass
        
    # if the method is UPDATE, update the data in the table with the respective columns
    elif method == "UPDATE":
        try:
            # update the data in the table where the primary_id is the same as the entry_id
            update_entry = '''UPDATE {} 
                            SET source = %s, primary_id = %s, gene_name = %s, organism = %s, expression_system = %s,
                            mitochondrial = %s, protein_sequence = %s, nucleotide_id = %s, nucleotide_sequence = %s, 
                            WHERE primary_id = %s'''.format(table)
            entry = (source, entry_id, gene_name, organism, expression_system, mitochondrial, protein_sequence, nucleotide_id, nucleotide_sequence)
            
            # execute the query and commit the changes to the database
            cursor.execute(update_entry, entry)
            DB.commit()
        except:
            logging.error(f'Entry {entry_id} not successfully updated in {table}!')
            pass
    
def insert_pandas(df, source, entry_id, gene_name, organism, expression_system, mitochondrial, protein_sequence, nucleotide_id, nucleotide_sequence):
    """Inserts the data into a pandas DataFrame."""
    try:
        # append the data to the pandas DataFrame with the respective columns
        # ignore the index to avoid the index column in the csv file
        # the dataframe can be empty or already contain data, but it has to have the same columns as the data to be inserted
        df = df._append({"source": source, "primary_id": entry_id, "gene_name": gene_name, "organism": organism, "expression_system": expression_system, "mitochondrial": mitochondrial, 
                        "protein_sequence": protein_sequence, "nucleotide_id": nucleotide_id, "nucleotide_sequence": nucleotide_sequence}, ignore_index=True)
        
    except:
        logging.error(f'Entry {entry_id} not successfully inserted in pandas DataFrame!')
        pass
    
    return df

def retrieve_sequence(ncbiDB="nuccore", entryID=None, rettype="fasta", retmode="text"):
    """Retrieves the nucleotide sequence from the NCBI database."""
    
    # since the script is run in parallel, ncbi only allows a certain number of requests per second
    # by providing an email address, 10 requests per second are allowed (so the script can run parallel on 10 cores)
    Entrez.email = mail
    Entrez.api_key = api_key
    
    # fetch the data from the NCBI database with the respective entryID
    handle = Entrez.efetch(db=ncbiDB, id=entryID, rettype=rettype, retmode=retmode)
    
    return handle.read()

def get_base_data(uniprotID):
    """Get the base data for the protein from the uniprot database."""
    try:
        # set up rest api query to get the data from the uniprot database
        query_url = 'https://rest.uniprot.org/uniprotkb/search?query=%s&fields=gene_primary,sequence,organism_name,cc_subcellular_location' %uniprotID
        
        # get the response from the rest api query
        response = session.get(url=query_url)
        
        # if the response is successful, parse the data from the response
        if response.status_code == 200:
            response = json.loads(response.text)
            
            try:
                # get the oganism if available for the uniprotID
                organism = response['results'][0]['organism']['scientificName']
            except:
                organism = "NaN"
                
            try:
                # get the gene name if available for the uniprotID
                gene_name = response['results'][0]['genes'][0]['geneName']['value']
            except:
                gene_name = "NaN"
            
            try:
                # get the sequence if available for the uniprotID
                sequence = response['results'][0]['sequence']['value']
            except: 
                sequence = "NaN"
    except:
        logging.error(f"Base data for protein {uniprotID} not available!")
        pass
    
    return organism,gene_name,sequence

def filter_sequence(sequences, searchID):
    """Filter the sequences based on the search criteria."""
    # the cds from ncbi includes multiple gene sequences, so the right nucleotide sequence 
    # has to be filtered containing the same cds_id (gene id from embl or refseq) as the protein sequence
    
    # set the search criteria to filter the sequences
    # create a set to store the unique sequences, and a list to store the matching data
    search_criteria = searchID
    unique_sequences = set()
    matching_data = []
    entries = sequences.split('>') # split the sequences

    for entry in entries:
        # loop through the entries and check if the search criteria is in the entry
        # it's possible to check for multiple criteria, that can be present in the entry or header 
        if any(criterion in entry for criterion in search_criteria):
            # get the header and the sequence from the entry
            header, sequence = entry.split('\n', 1)
            sequence = sequence.strip()
            
            # Check if the sequence is already added to ensure uniqueness
            # this results in having only one sequence per cds_id with the correct nucleotides / codons
            if sequence not in unique_sequences:
                unique_sequences.add(sequence)
                matching_data.append((header, sequence.replace('\n',"")))

    return matching_data

def check_isoform(uniprotID):
    """Check if the protein in uniprot has available isoforms."""
    # for some eukaryotic entries, uniprot provides isoforms that derive from for instance alternative splicing
    # this function checks if the protein has isoforms and returns True if it has, otherwise False
    
    # set up the rest api query to get the isoform data (alternative products) from the uniprot database
    url_iso = 'https://rest.uniprot.org/uniprotkb/search?query=%s&fields=cc_alternative_products' %uniprotID

    # get the response from the rest api query with the isoforms
    response_gene = session.get(url_iso)

    if response_gene.status_code == 200: 
        response = json.loads(response_gene.text)
        
        try:
            # check if the isoform data is available for the protein
            if any(response['results'][0]['comments'][0]['isoforms']):
                return True
        except:
            return False

def get_isoform(uniprotID):
    """Get the isoform data for the protein from the uniprot database."""
    # if there's isoform data available for the protein, the isoform data is retrieved
    # sometimes multiple isoforms are available for a protein, which also will be retrieved
    
    # set up the rest api query to get the isoform data from the uniprot database
    # get the cross references for the isoforms from the uniprot database
    url_xref = 'https://rest.uniprot.org/uniprotkb/search?query=%s&fields=xref_refseq' %uniprotID

    # get the response from the rest api query with the isoform data
    response = session.get(url=url_xref)

    if response.status_code == 200:
        response = json.loads(response.text)
        
        # if the response is successful, parse the data from the response
        if response['results'][0]['uniProtKBCrossReferences']:
            # get the protein ids, cds ids, and isoform ids for the isoforms and store them in lists
            protein_ids = []
            cds_ids = []
            isoformIds = []
            
            for cds_ref in response['results'][0]['uniProtKBCrossReferences']:
                # loop through the cross references and get the protein id, cds id, and isoform id
                protein_id = cds_ref['id']
                
                # check if the protein id is a refseq protein id annotated with NP
                if protein_id.startswith("NP"):
                    protein_ids.append(protein_id)
                    
                    # the the corresponding cds_id (refseq) for the protein id
                    cds_id = cds_ref['properties'][0]['value']

                    # check if the cds id is a refseq id annotated with NM
                    if cds_id.startswith("NM"):
                        cds_ids.append(cds_id)
                    
                try:
                    # also get the isoform id if available
                    isoformId = cds_ref['isoformId']
                    isoformIds.append(isoformId)
                except:
                    logging.error(f"Isoform ID for protein {uniprotID} not available!")
                    pass
                    
            return protein_ids, cds_ids, isoformIds

def get_cds(uniprotID):
    """Get the nucleotide coding sequences for the protein from the uniprot database."""
    # set up rest api query to get the cds data from the uniprot database
    # fetching the embl cross reference for non-isoform proteins is easier
    # but is equivalent to the refseq cross references
    # this functions is thus used for non-isoform proteins / entries
    url_cds = 'https://rest.uniprot.org/uniprotkb/search?query=%s&fields=xref_embl' %uniprotID

    # get the response from the rest api query with the cds data
    response_cds = session.get(url_cds)

    if response_cds.status_code == 200: 
        # if the response is successful, parse the data from the response
        # get the cds ids and protein ids for the cds and store them in lists
        cds_ids = []
        protein_ids = []
        
        response = json.loads(response_cds.text)
        
        # check if the response contains cross references for the cds
        if response['results'][0]['uniProtKBCrossReferences']:
            for cds_ref in response['results'][0]['uniProtKBCrossReferences']:
                # loop over the cross references and get the cds id and protein id
                cds_id = cds_ref['id']
                
                for properties in cds_ref['properties']:
                    # the cds reference contain properties and can be filtered by the values of their keys
                    
                    # get the protein id and add it to the list if the key is ProteinId
                    # add them only if their are no sequence conflicts like "sequence problems", "different initiation", etc.
                    if properties['key'] == 'ProteinId':
                        protein_id = properties['value']
                
                    if properties['key'] == 'Status':
                        status_code = properties['value']
                
                # only add the cds id and protein id if the status code is not "-"
                if status_code == "-":
                    cds_ids.append(cds_id)
                    protein_ids.append(protein_id)
                        
        return cds_ids, protein_ids
    
def main():
    parser = argparse.ArgumentParser(
        prog="UPmapNT",
        description="Map uniprot ID to nucleotide sequence and prints an allignment of the pdb protein sequence to the nucleotide sequence"
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
        '--create-fasta', choices=['protein', 'nucleotide', 'all'], dest="create_fasta", default=False,
        help="Create a fasta file with the nucleotide sequences, protein sequences or both."
    )
    parser.add_argument(
        '-w', '--overwrite', action="store_true", dest="overwrite", required=False, default=False,
        help='If file name already exists, overwrite it. Default is False.' 
    )
    parser.add_argument(
        '-m', '--mail', type=str, dest="api_mail", required=True, default=False,
        help='Provide mail from ncbi api account: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/' 
    )
    parser.add_argument(
        '-k', '--key', type=str, dest="api_key", required=True, default=False,
        help='Provide api key from ncbi api account: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/' 
    )
    args = parser.parse_args()
    
    # connect to the SQL database if the respective flag is set /* might be removed or extended in future */
    if args.sql:
        nustruDB = connect_DB()
        
    # alternatively, store the data in a pandas DataFrame amnd create a csv file with the column names to store the data
    elif args.pandas:
        nucleotide_protein_seqs_df = pd.DataFrame(columns=["source", "primary_id", "gene_name", "organism", "expression_system", "mitochondrial", "protein_sequence", "nucleotide_id", "nucleotide_sequence"])
        if Path(f'{args.output_path}/{args.name}.csv').exists() or args.overwrite:
                nucleotide_protein_seqs_df.to_csv(f'{args.output_path}/{args.name}.csv', mode='w', index=False, header=True)
        else:
            print(f"Error: {args.output_path} already exists. Use -w to overwrite.")
            exit(1)
            
    else:
        print("Error: No output format specified. Use --sql or --pandas.")
        exit(1) 
        
    logging.basicConfig(filename=f'{args.output_path}/{args.name}.log',
            filemode='a',
            format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
            datefmt='%H:%M:%S',
            level=logging.ERROR)
    
    # set the api informations as global variables to use them in the functions
    # providing faster access to the api, so should be given
    global mail, api_key
    mail = args.api_mail
    api_key = args.api_key
    
    with open(args.entryID,'r') as entryIDs_file:
        # read the file with the uniprot IDs and split them by comma
        entryIDs_file = entryIDs_file.read()
        
        # create a progress bar to show the progress of the mapping
        progress = tqdm(total=len(entryIDs_file.replace("\n", ",").split(',')))

        for uniprotID in entryIDs_file.replace("\n", ",").split(','):
            # loop through the uniprot IDs, separate them by comma, and get the necessary data
            
            # assign the base data to the respective variables
            organism, gene_name, sequence = get_base_data(uniprotID)
    
            try:
                try:
                    # check if the protein has isoforms and get the isoform data if available
                    if check_isoform(uniprotID):
                        # get the isoform data for the protein
                        protein_ids, cds_ids, isoform_ids = get_isoform(uniprotID)
                        
                        # store the isoform protein sequences in a list
                        isoform_protein_sequences = []
                        for cds_id, protein_id, isoform_id in zip(cds_ids, protein_ids, isoform_ids):
                            # loop through the isoform data and get the nucleotide sequence for the isoform
                            # the cds nucleotide sequence can be directly retrieved by using "fasta_cds_na" as rettype
                            # find the right nucleotide sequence for the protein sequence with the same cds_id
                            nt_response_sequence = retrieve_sequence(entryID=cds_id, rettype="fasta_cds_na")
                            nt_response_sequence = filter_sequence(sequences=nt_response_sequence,searchID=[cds_id])
                            
                            # get the protein sequence for the isoform from ncbi
                            protein_response_sequence = retrieve_sequence(ncbiDB="protein",entryID=protein_id, rettype="fasta")
                            protein_response_sequence = filter_sequence(sequences=protein_response_sequence,searchID=[protein_id])
                            isoform_protein_sequences.append(protein_response_sequence[0][1])
                            
                            # create fasta files with the nucleotide and / or protein sequences of the PDB entries if the respective flags are set
                            if args.create_fasta == "nucleotide" or args.create_fasta == "all":
                                with open(f'{args.output_path}/{args.name}.fasta', 'a') as f:
                                    f.write(f'>{isoform_id} | {cds_id} {organism}\n{nt_response_sequence[0][1]}\n')
                            if args.create_fasta == "protein" or args.create_fasta == "all":
                                with open(f'{args.output_path}/{args.name}.fasta', 'a') as f:
                                    f.write(f'>{isoform_id} | {protein_id} {organism}\n{protein_response_sequence[0][1]}\n')
                    
                            # store the data in the SQL database /* might be removed or extended in future */
                            # or pandas DataFrame if the respective flags are set
                            if args.sql:
                                execute_database(DB=nustruDB, method="INSERT", table="nucleotide_protein_seqs", source="uniprot", entry_id=isoform_id, gene_name=gene_name, organism=organism, 
                                                expression_system="NaN", mitochondrial="False", protein_sequence=sequence,
                                                nucleotide_id=cds_id, nucleotide_sequence=nt_response_sequence[0][1])
                            else:
                                nucleotide_protein_seqs_df = insert_pandas(df=nucleotide_protein_seqs_df, source="uniprot", entry_id=isoform_id, gene_name=gene_name, organism=organism, 
                                                expression_system="NaN", mitochondrial="False", protein_sequence=sequence,
                                                nucleotide_id=cds_id, nucleotide_sequence=nt_response_sequence[0][1])
                                
                                # prevent loss of data by writing the data sequentially to the csv file
                                # delete the pandas DataFrame to free up memory
                                nucleotide_protein_seqs_df.to_csv(f'{args.output_path}/{args.name}.csv', mode='a', index=False, header=False)
                                del nucleotide_protein_seqs_df
                                
                except:
                    logging.error(f"Isoform data and nucleotide sequence for protein {uniprotID} not available!")
                    continue
                 
                try:        
                    # if the protein has no isoforms, get the cds data for the protein
                    if not check_isoform(uniprotID) or not sequence in isoform_protein_sequences:
                        # get the cds data for the protein
                        cds_ids, protein_ids = get_cds(uniprotID)
                        
                        for cds_id, protein_id in zip(cds_ids, protein_ids):
                            # loop through the cds data and get the nucleotide sequence for the protein
                            # check if the retrieved protein sequence is the same as the protein sequence
                            if retrieve_sequence(ncbiDB="protein", entryID=protein_id).partition("\n")[2].strip().replace('\n','') == sequence:
                                # for same protein sequences, get the nucleotide sequence for the protein
                                response_sequences = retrieve_sequence(entryID=cds_id, rettype="fasta_cds_na")
                                
                                # filter the nucleotide sequence based on the search criteria, the header should contain the uniirotID or the protein_id
                                matched_id, matched_seq = filter_sequence(sequences=response_sequences,searchID=[uniprotID, protein_id])[0]
                                
                                # create fasta files with the nucleotide and / or protein sequences of the PDB entries if the respective flags are set
                                if args.create_fasta == "nucleotide" or args.create_fasta == "all":
                                    with open(f'{args.output_path}/{args.name}.fasta', 'a') as f:
                                        f.write(f'>{uniprotID} | {cds_id} {organism}\n{nt_response_sequence[0][1]}\n')
                                if args.create_fasta == "protein" or args.create_fasta == "all":
                                    with open(f'{args.output_path}/{args.name}.fasta', 'a') as f:
                                        f.write(f'>{uniprotID} | {protein_id} {organism}\n{protein_response_sequence[0][1]}\n')
                        
                                # store the data in the SQL database /* might be removed or extended in future */
                                # or pandas DataFrame if the respective flags are set
                                if args.sql:
                                    execute_database(DB=nustruDB, method="INSERT", table="nucleotide_protein_seqs", source="uniprot", entry_id=uniprotID, gene_name=gene_name, organism=organism, 
                                                    expression_system="NaN", mitochondrial="False", protein_sequence=sequence,
                                                    nucleotide_id=cds_id, nucleotide_sequence=matched_seq)
                                else:
                                    nucleotide_protein_seqs_df = insert_pandas(df=nucleotide_protein_seqs_df, source="uniprot", entry_id=uniprotID, gene_name=gene_name, organism=organism, 
                                                    expression_system="NaN", mitochondrial="False", protein_sequence=sequence,
                                                    nucleotide_id=cds_id, nucleotide_sequence=matched_seq)
                                    
                                    # prevent loss of data by writing the data sequentially to the csv file
                                    # delete the pandas DataFrame to free up memory
                                    nucleotide_protein_seqs_df.to_csv(f'{args.output_path}/{args.name}.csv', mode='a', index=False, header=False)
                                    del nucleotide_protein_seqs_df
                                
                                # break the loop if a matching protein sequence is found (with no sequence conflicts)
                                break
                except:
                    logging.error(f"Nucleotide sequence for protein {uniprotID} not available!")
                    continue

            except:
                logging.error(f"No isoform or nucleotide sequence for protein {uniprotID} available!")
                continue
            
            # create a new pandas DataFrame to store the data
            nucleotide_protein_seqs_df = pd.DataFrame(columns=["source", "primary_id", "gene_name", "organism", "expression_system", "mitochondrial", "protein_sequence", "nucleotide_id", "nucleotide_sequence"])
            # update the progress bar
            progress.update(1)
            
if __name__ == '__main__':
    main()