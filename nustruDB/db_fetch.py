#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
from biopandas.pdb import PandasPdb 
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP


df = pd.read_csv('/Users/dominiquefastus/master_project/NuStru/nustruDB/E_COLI_K12_02_reduced.csv', nrows=10, index_col=False)

def fetch_pdb_and_plddt(data):
    if data['source'] == 'uniprot':
        ppdb = PandasPdb().fetch_pdb(uniprot_id=data['primary_id'], source="alphafold2-v2")
        residue_b_factors = ppdb.df["ATOM"][["b_factor","residue_number"]]


        residue_b_factors = residue_b_factors.groupby(["residue_number"]).mean().round(4)

        residue_b_factors_dict = residue_b_factors.to_dict()
        
        
        ppdb.to_pdb(path='/Users/dominiquefastus/Downloads/{}.pdb'.format(data['primary_id']))
    
df.apply(fetch_pdb_and_plddt, axis=1)
        

"""
pdb_ids = df.loc[df['source'] == 'pdb']
uniprot_ids = df.loc[df['source'] == 'uniprot']

pdb_ids = pdb_ids['primary_id'].tolist()
uniprot_ids = uniprot_ids['primary_id'].tolist()

url = '//www.rcsb.org/pdb/files/fasta.txt?structureIdList=%s' % ','.join(pdb_ids)

response = requests.get(url=url)

if response.status_code == 200: 
    with open
"""