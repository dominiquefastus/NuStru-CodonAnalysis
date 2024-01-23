import argparse
import requests
import json
import pandas as pd

from bs4 import BeautifulSoup

from dataclasses import dataclass
from string import Template

from Bio import Entrez

def get_allignment(entryID, db="pdb"):
    allignment_range = {}
    
    if db == "pdb":
        pdbID = entryID
        
    elif len(entryID) > 5:
        uniprotID = entryID
        
    else:
        raise Exception()
        print("Not a valid id")
        exit(1)
        
    # saf

    query = '{alignment(from: PDB_ENTITY, to: NCBI_GENOME, queryId: "%s" ) { query_sequence target_alignment { target_id orientation aligned_regions { query_begin query_end target_begin target_end } } }}' % pdbID
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
            range = [position['target_begin'], position['target_end']]
            allignment_range[id] = range
      
    return pdb_sequence, genomeID, oritentation, allignment_range
    
        
def map_uniprot(uniprotID):
    pass
    
def retrieve_nucleotide_seqs(genomeID=None, seqSTART=None, seqEND=None, orientation=None):
    # same as requesting by url
    # f"https://www.ncbi.nlm.nih.gov/nuccore/{genomeID}?report=fasta&from={seqSTART}&to={seqEND}"
    
    # set strand orientation (positive or negative strand)
    if orientation == 1:
        orientation = 1
    else:
        orientation = 2
        
    Entrez.email = "dominique.fastus@biochemistry.lu.se"    
    handle = Entrez.efetch(db="nuccore", id=genomeID, seq_start=seqSTART, seq_stop=seqEND, strand=orientation, rettype="fasta", retmode="text")
    
    return handle.read().partition("\n")[2].strip()
    
def console(pdb_seq, nt_seq):
    pass

def main():
    parser = argparse.ArgumentParser(
        prog="PDBmapNT",
        description="Maps PDB ID to nucleotide sequence and prints an allignment of the pdb protein sequence to the nucleotide sequence"
    )
    parser.add_argument('entryID')
    args = parser.parse_args()
    
    pdb_sequence, genomeID, oritentation, allignment_range = get_allignment(args.entryID)
    
    print(pdb_sequence, genomeID, oritentation, allignment_range)
    
    nu_sequence = ""
    for seqSTART, seqEND in allignment_range.values():
        nu_sequence += retrieve_nucleotide_seqs(genomeID=genomeID, seqSTART=seqSTART, seqEND=seqEND)
    
    print(nu_sequence.replace('\n',''))
    
    

if __name__ == '__main__':
    main()