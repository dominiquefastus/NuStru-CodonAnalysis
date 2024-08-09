#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: PDBmapNT [-h] -i ENTRYID [--sql] [--pandas] -o OUTPUT_PATH -n NAME [--map-uniprot] [--create-fasta {protein,nucleotide,all}] [-w]

Map PDB ID to nucleotide sequence and prints an allignment of the pdb protein sequence to the nucleotide sequence

options:
  -h, --help            show this help message and exit
  -i ENTRYID, --input ENTRYID
                        Path to file containing PDB IDs with comma separation.
  --sql                 Store the data in a SQL database.
  --pandas              Store the data in a pandas DataFrame.
  -o OUTPUT_PATH, --output OUTPUT_PATH
                        Output path or directory to store the log file and the data.
  -n NAME, --name NAME  Name of the output files and log file.
  --map-uniprot         Map uniprot ID to nucleotide sequence.
  --create-fasta {protein,nucleotide,all}
                        Create a fasta file with the nucleotide sequences, protein sequences or both.
  -w, --overwrite       If file name already exists, overwrite it. Default is False.
  
The script will create a csv file where the data from the PDB ID mapping is stored. It is possible to also map the respective uniprot entries
that are associated with the PDB ID. The script will also create fasta files with the nucleotide and protein sequences of the PDB and uniprot if the respective flags are set.
Example: python PDBmapNT.py -i Example/examples_nustruDB/example_pdbIDs.txt --pandas -o . -n example_pdbIDs_nustru --create-fasta all -w
"""

import argparse
import requests
import logging
import json

import pandas as pd

from tqdm import tqdm
from pathlib import Path

import mysql.connector
from getpass import getpass
from requests.adapters import HTTPAdapter, Retry

# create a session with retries, with a total of 6 retries, a backoff factor of 0.2, and status codes 502, 503, and 504
session = requests.Session()
retries = Retry(total=6, backoff_factor=0.2, status_forcelist=[ 502, 503, 504 ])
session.mount('https://', HTTPAdapter(max_retries=retries)) # mount the session with the retries

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
       logging.error("Database connection is not established.")
       exit(1)

    cursor = DB.cursor() # create a cursor object to execute SQL queries

    # if the method is INSERT, insert the data into the table with the respective columns
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

def get_base_data(entryID):
    """Get the base data for the protein from the PDB database."""
    # by using the polymer entity of the entry id, it is possible to get the data for each individual chain of the protein
    # this works also for monomers as then only the data for one chain will be fetched
    try:
        # set up graphql query to get the organism, gene name, and expression system of the protein
        # the query is based on the entry ID of the protein and retrieves the entity ids (chains) of the protein
        #  query = '{entry(entry_id: "%s")  { rcsb_id polymer_entities { rcsb_entity_source_organism { scientific_name ncbi_scientific_name rcsb_gene_name { value } } polymer_entity_instances { rcsb_polymer_instance_feature { type feature_id feature_positions { end_seq_id beg_seq_id } } } rcsb_entity_host_organism { ncbi_scientific_name } } rcsb_entry_container_identifiers { polymer_entity_ids } } }' % entryID
        query = '{entry(entry_id: "%s") { rcsb_id polymer_entities { rcsb_entity_source_organism { scientific_name ncbi_scientific_name rcsb_gene_name { value } } rcsb_entity_host_organism { ncbi_scientific_name } } rcsb_entry_container_identifiers { polymer_entity_ids } }}' % entryID
        url = f'https://data.rcsb.org/graphql?query={query}'
        
        # get the response from the graphql query
        response = session.get(url=url)

        # if the query or response is successful, parse the data from the fetched json
        if response.status_code == 200:
            response = json.loads(response.text)

            try:
                # get the organism if available for the entity
                organism = response['data']['entry']['polymer_entities'][0]['rcsb_entity_source_organism'][0]['ncbi_scientific_name']
            except:
                organism = "NaN"
        
            try:
                # get the gene name if available for the entity
                gene_name = response['data']['entry']['polymer_entities'][0]['rcsb_entity_source_organism'][0]['rcsb_gene_name'][0]['value']
            except:
                gene_name = "NaN"
                
            try:
                # get the expression system if available for the entity
                expression_system = response['data']['entry']['polymer_entities'][0]['rcsb_entity_host_organism'][0]['ncbi_scientific_name']
            except:
                expression_system = "NaN"
                
            try:
                # get the polymer entity ids of the protein (chains)
                # create a list of the polymer entity ids with the entry id as prefix (chain ids)
                polymer_entity_ids = response['data']['entry']['rcsb_entry_container_identifiers']['polymer_entity_ids']
                pdb_entity_id_list = [entryID + '_' + entity_id for entity_id in polymer_entity_ids]
            except: 
                # if the polymer entity ids are not available, set the variables to None
                polymer_entity_ids = None
                pdb_entity_id_list = None
    except:
        logging.error(f"Base data for protein {entryID} not available!")
        pass
    
    return organism, gene_name, expression_system, polymer_entity_ids, pdb_entity_id_list

def map_entity_to_chains(entites, entryID):
    """Map the polymer entity to the chains of the protein."""    
    # the polymer entity is mapped to the chains of the protein by using the pdbx_strand_id
    pdb_chains_entity_list = []
    for entity in entites:
        # loop over the entities and get the pdbx_strand_id of the entity   
        query = '{polymer_entity(entity_id: "%s", entry_id: "%s") { rcsb_id entity_poly { pdbx_strand_id } }}' % (entity, entryID)
        url = f'https://data.rcsb.org/graphql?query={query}'
        
        response = session.get(url=url)
        
        if response.status_code == 200:
            # parse the response and get the pdbx_strand_id of the entity
            response = json.loads(response.text)
            pdbx_strand_id = response['data']['polymer_entity']['entity_poly']['pdbx_strand_id']
            
            # create a list of the chains of the protein with the entry id as prefix
            pdb_chain_entity = entryID + '_' + pdbx_strand_id.split(',')[0]
            pdb_chains_entity_list.append(pdb_chain_entity)
    
    return pdb_chains_entity_list
       
def get_allignment(entryID, id_type):
    """Get the allignment of the protein sequence to the nucleotide sequence."""
    try:
        allignment_range = {}
        exon_shift_range = {}
        
        if id_type == "pdb":
            # get the allignment of the protein sequence to the nucleotide sequence by the pdb entity id
            pdbID = entryID
            query = '{alignment(from: PDB_ENTITY, to: NCBI_GENOME, queryId: "%s" ) { query_sequence target_alignment { target_id orientation aligned_regions { query_begin query_end target_begin target_end exon_shift } } }}' % pdbID
            
        elif id_type == "uniprot":
            # if the corresponding uniprot id wants to be fetched, the uniprot id is used to get the allignment
            uniprotID = entryID
            query = '{alignment(from: UNIPROT, to: NCBI_GENOME, queryId: "%s" ) { query_sequence target_alignment { target_id orientation aligned_regions { query_begin query_end target_begin target_end exon_shift } } }}' % uniprotID
            
        else:
            logging.error(f"ID type {id_type} not recognized!")
            pass
            
        # define graphql url from pdb
        url = f'https://1d-coordinates.rcsb.org/graphql?query={query}'
        
        response = session.get(url=url)
        
        if response.status_code == 200:
            # make the graphql query and parse the response
            response = json.loads(response.text)
            
            # get the protein sequence, genome ID, orientation, allignment range, and exon shift range
            # these information are needed to retrieve the nucleotide sequence
            pdb_sequence = response['data']['alignment']['query_sequence']
            genomeID = response['data']['alignment']['target_alignment'][0]['target_id']
            oritentation = response['data']['alignment']['target_alignment'][0]['orientation']
            
            # loop through the sequence positions
            for id, position in enumerate(response['data']['alignment']['target_alignment'][0]['aligned_regions']):
                # for each sequence position, get the target begin and end positions and the exon shift
                range = [position['target_begin'], position['target_end']]
                exon_shift_range[id] = position['exon_shift']
                allignment_range[id] = range
    except:
        # if the allignment of the protein sequence to the nucleotide sequence is not available, skip the entry
        pass
      
    return pdb_sequence, genomeID, oritentation, allignment_range, exon_shift_range
    
        
def map_uniprot(pdbID):
    """Map the pdb ID to the uniprot ID."""
    # sometime the uniprot id linked to the pdb entry wants to be fetched
    # align the pdb entity to the uniprot entity to get the uniprot id
    query = '{alignment(from: PDB_ENTITY, to: UNIPROT, queryId: "%s" ) { query_sequence target_alignment { target_id target_sequence } }}' % pdbID
    # define graphql url from pdb
    url = f'https://1d-coordinates.rcsb.org/graphql?query={query}'
    
    response = session.get(url=url)
    
    if response.status_code == 200:
        # extract the uniprot id from the response
        response = json.loads(response.text)
        uniprot_id = response['data']['alignment']['target_alignment'][0]['target_id']
        
    return uniprot_id
    
def retrieve_nucleotide_seqs(genomeID=None, orientation=None, seqSTART=None, seqEND=None):
    """Retrieve the nucleotide sequence from the genome.""" 
    # same as requesting by url
    # f"https://www.ncbi.nlm.nih.gov/nuccore/{genomeID}?report=fasta&from={seqSTART}&to={seqEND}"
    
    # check if the start and end positions are in the right order
    # if not, switch the positions
    if seqSTART > seqEND:
        seqSTART, seqEND = seqEND, seqSTART
        
    # set strand orientation (positive or negative strand)
    # entrez requires 1 for positive and 2 for negative strand
    if orientation == -1:
        orientation = 2
    else:
        orientation = 1
           
    # fetch the nucleotide sequence from the genome by the genome ID, orientation, start and end positions
    # return the sequence in fasta format and text mode, also remove the first line of the fasta file and strip the sequence
    handle = Entrez.efetch(db="nuccore", id=genomeID, seq_start=seqSTART, seq_stop=seqEND, strand=orientation, rettype="fasta", retmode="text")
    
    return handle.read().partition("\n")[2].strip()

def main():
    parser = argparse.ArgumentParser(
        prog="PDBmapNT",
        description="Map PDB ID to nucleotide sequence and prints an allignment of the pdb protein sequence to the nucleotide sequence"
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
        '--map-uniprot', action="store_true", dest="map_uniprot", default=False,
        help="Map uniprot ID to nucleotide sequence."
    )
    parser.add_argument(
        '--create-fasta', choices=['protein', 'nucleotide', 'all'], dest="create_fasta", default=False,
        help="Create a fasta file with the nucleotide sequences, protein sequences or both."
    )
    parser.add_argument(
        '-w', '--overwrite', action="store_true", dest="overwrite", required=False, default=False,
        help='If file name already exists, overwrite it. Default is False.' 
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
        
    with open(args.entryID,'r') as entryIDs_file:
        # read the file with the PDB IDs and split the IDs by comma
        entryIDs_file = entryIDs_file.read()
        
        # create a progress bar to show the progress of the mapping
        progress = tqdm(total=len(entryIDs_file.replace(" ", "").split(',')))

        # loop over the PDB IDs and get the base data for the protein
        for pdb_id in entryIDs_file.replace(" ", "").split(','):
            # assign base data to variables, that will be stored in the database or pandas DataFrame
            # get the chains from the polymer entity ids, including multiple chains if available
            organism, gene_name, expression_system, polymer_entity_ids, pdb_entity_id_list = get_base_data(entryID=pdb_id)
            pdb_chains_entity_list = map_entity_to_chains(entites=polymer_entity_ids, entryID=pdb_id)
            
            # loop over the proteins and its chains and get the allignment of the protein sequence to the nucleotide sequence
            for pdb_id, pdb_chain_id in zip(pdb_entity_id_list, pdb_chains_entity_list):
                pdb_entry = pdb_chain_id

                # most pdb entries have mapped uniprot entries, which can be fetched
                # if the respective flag is set, the uniprot id will be fetched
                if args.map_uniprot:
                    try:
                        uniprotID = map_uniprot(pdb_id)
                    except:
                        logging.error(f'No uniprot ID for protein {pdb_id} available')
                        uniprotID = None
                    
                try:
                    # get the information from the allignment of the protein sequence to the nucleotide sequence
                    # this includes the query sequence with the aligned regions and the exon shift (if present)
                    pdb_sequence, genomeID, oritentation, allignment_range, exon_shift_range = get_allignment(entryID=pdb_id, id_type="pdb")
                    
                    # add the nucleotide sequence parts or ranges to one whole sequence
                    nu_sequence = ""
                    # loop over the allignment range and exon shift range to get the nucleotide sequence
                    # the nucleotide sequence is fetched by the genome ID, orientation, and the start and end positions
                    for(seqSTART, seqEND), (exon_range) in zip(allignment_range.values(), exon_shift_range.values()):
                        nu_sequence += retrieve_nucleotide_seqs(genomeID=genomeID, orientation=oritentation, seqSTART=seqSTART, seqEND=seqEND)
                        
                        # if an exon range is present, get the nucleotide sequence of the exon range
                        if len(exon_range) == 2:
                            exonSTART = exon_range[0]
                            exonEND = exon_range[1]
                            
                            nu_sequence += retrieve_nucleotide_seqs(genomeID=genomeID, orientation=oritentation, seqSTART=exonSTART, seqEND=exonEND)
                            
                        # if only one exon range is present, get the nucleotide sequence of the exon range
                        elif len(exon_range) == 1:
                            exonSTART = exon_range[0]
                            exonEND = exon_range[0]
                        
                            
                            nu_sequence += retrieve_nucleotide_seqs(genomeID=genomeID, orientation=oritentation, seqSTART=exonSTART, seqEND=exonEND)
                            
                        # if no exon range is present, skip the exon range
                        elif len(exon_range) == 0:
                            pass
                        
                        else:
                            print(exon_range)
                        
                    # remove the newline characters from the nucleotide sequence
                    # make a complete sequence of the nucleotides / codons
                    nu_sequence = nu_sequence.replace('\n','')
                    
                    # create fasta files with the nucleotide and / or protein sequences of the PDB entries if the respective flags are set
                    if args.create_fasta == "nucleotide" or args.create_fasta == "all":
                        with open(f'{args.output_path}/{args.name}_nt_pdb.fasta', 'a') as f:
                            f.write(f'>{pdb_entry}| {genomeID} [{allignment_range}] [{exon_shift_range}] {organism}\n{nu_sequence}\n')
                    if args.create_fasta == "protein" or args.create_fasta == "all":
                        with open(f'{args.output_path}/{args.name}_prot_pdb.fasta', 'a') as f:
                            f.write(f'>{pdb_entry}| {organism}\n{pdb_sequence}\n')
                        
                    # store the data in the SQL database /* might be removed or extended in future */
                    # or pandas DataFrame if the respective flags are set
                    if args.sql:
                        execute_database(DB=nustruDB, method="INSERT", table="nucleotide_protein_seqs", source="pdb", entry_id=pdb_entry, gene_name=gene_name, organism=organism, 
                                        expression_system=expression_system, mitochondrial="False", protein_sequence=pdb_sequence,
                                        nucleotide_id=genomeID, nucleotide_sequence=nu_sequence)
                    else:
                        nucleotide_protein_seqs_df = insert_pandas(df=nucleotide_protein_seqs_df, source="pdb", entry_id=pdb_entry, gene_name=gene_name, organism=organism, 
                                                                expression_system=expression_system, mitochondrial="False", protein_sequence=pdb_sequence,
                                                                nucleotide_id=genomeID, nucleotide_sequence=nu_sequence)
                        
                        # prevent loss of data by writing the data sequentially to the csv file
                        # delete the pandas DataFrame to free up memory
                        nucleotide_protein_seqs_df.to_csv(f'{args.output_path}/{args.name}.csv', mode='a', index=False, header=False)
                        del nucleotide_protein_seqs_df
                        
                except:
                    logging.error(f"Genomic coordinates for protein {pdb_id} not available!")
                    continue
                
                # create a new pandas DataFrame to store the data
                nucleotide_protein_seqs_df = pd.DataFrame(columns=["source", "primary_id", "gene_name", "organism", "expression_system", "mitochondrial", "protein_sequence", "nucleotide_id", "nucleotide_sequence"])
                    
                # map the uniprot ID from the corresponding PDB entries to the nucleotide sequence if the respective flag is set
                if args.map_uniprot:
                    try:
                        # get the information from the allignment of the protein sequence to the nucleotide sequence
                        uniprot_sequence, genomeID, oritentation, allignment_range, exon_shift_range = get_allignment(entryID=uniprotID, id_type="uniprot")
                        
                        # identical as before loop over the allignment range and exon shift range to get the nucleotide sequence
                        # uniprot entries probably won't have exon ranges, but the loop is kept for consistency
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
                        
                        # remove the newline characters from the nucleotide sequence
                        nu_sequence = nu_sequence.replace('\n','')
                        
                        # create fasta files with the nucleotide and / or protein sequences of the uniprot entries if the respective flags are set
                        if args.create_fasta == "nucleotide" or args.create_fasta == "all":
                            with open(f'{args.output_path}/{args.name}_nt_uniprot.fasta', 'a') as file:
                                file.write(f'>{uniprotID}| {genomeID} [{allignment_range}] [{exon_shift_range}] {organism}\n{nu_sequence}\n')
                        if args.create_fasta == "protein" or args.create_fasta == "all":
                            with open(f'{args.output_path}/{args.name}_prot_uniprot.fasta', 'a') as file:
                                file.write(f'>{uniprotID}| {organism}\n{uniprot_sequence}\n')
                        
                        # store the data in the SQL database /* might be removed or extended in future */
                        # or pandas DataFrame if the respective flags are set
                        if args.sql:
                            execute_database(DB=nustruDB, method="INSERT", table="nucleotide_protein_seqs", source="uniprot", entry_id=uniprotID, gene_name=gene_name, organism=organism, 
                                            expression_system=expression_system, mitochondrial="False", protein_sequence=uniprot_sequence,
                                            nucleotide_id=genomeID, nucleotide_sequence=nu_sequence)
                        else:
                            nucleotide_protein_seqs_df = insert_pandas(df=nucleotide_protein_seqs_df, source="uniprot", entry_id=uniprotID, gene_name=gene_name, organism=organism, 
                                                                    expression_system=expression_system, mitochondrial="False", protein_sequence=uniprot_sequence,
                                                                    nucleotide_id=genomeID, nucleotide_sequence=nu_sequence)
                            
                            # store sequentially the data in the csv file and delete the pandas DataFrame to free up memory
                            nucleotide_protein_seqs_df.to_csv(f'{args.output_path}/{args.name}.csv', mode='a', index=False, header=False)
                            del nucleotide_protein_seqs_df
                            
                    except:
                        logging.error(f"Genomic coordinates for protein {uniprotID} not available!")
                        continue
                
                # create a new pandas DataFrame to store the data
                nucleotide_protein_seqs_df = pd.DataFrame(columns=["source", "primary_id", "gene_name", "organism", "expression_system", "mitochondrial", "protein_sequence", "nucleotide_id", "nucleotide_sequence"])
                # update the progress bar
                progress.update(1)            

if __name__ == '__main__':
    main()