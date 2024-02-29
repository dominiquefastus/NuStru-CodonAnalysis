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
        residue_plddt = ppdb.df["ATOM"][["b_factor","residue_number"]]


        residue_plddt = residue_plddt.groupby(["residue_number"]).mean().round(4)

        residue_plddt_dict = residue_plddt.to_dict()
        data['plddt'] = residue_plddt_dict
        
    else:
        ppdb = PandasPdb().fetch_pdb(pdb_code=data['primary_id'], source="pdb")
        residue_b_factors = ppdb.df["ATOM"][["b_factor","residue_number"]]
        
        residue_b_factors = residue_b_factors.groupby(["residue_number"]).mean().round(4)

        residue_b_factors_dict = residue_b_factors.to_dict()
        
        data['plddt'] = residue_b_factors_dict
    

    new_data = pd.DataFrame({'source': data['source'], 'primary_id': data['primary_id'], 'gene_name': data['gene_name'], 'organism': data['organism'], 'expression_system': data['expression_system'],
                            'protein_sequence': data['protein_sequence'], 'nucleotide_id': data['nucleotide_id'], 'nucleotide_sequence': data['nucleotide_sequence'], 'plddt': data['plddt']}, index=[0])
    new_data.to_csv('/Users/dominiquefastus/Downloads/structure_info.csv', mode='a', index=False, header=False)
    
    ppdb.to_pdb(path='/Users/dominiquefastus/Downloads/{}.pdb'.format(data['primary_id']))


with open('/Users/dominiquefastus/Downloads/structure_info.csv', mode='w') as f:
    f.write('source,primary_id,gene_name,organism,expression_system,protein_sequence,nucleotide_id,nucleotide_sequence,plddt,secondary_structure\n')
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