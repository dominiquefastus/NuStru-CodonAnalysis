#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq

df = pd.read_csv('/Users/dominiquefastus/master_project/NuStru/nustruDB/E_COLI_K12_02.csv', index_col=False)


df_reduced = df.drop_duplicates(subset=['primary_id','nucleotide_sequence'])

df_reduced.to_csv('/Users/dominiquefastus/master_project/NuStru/nustruDB/E_COLI_K12_02_reduced.csv', index=False)

def nucleotide_to_protein(data, codons = 3):
    translated_sequence = data['nucleotide_sequence'].str[:codons].apply(lambda x: Seq(x).translate())
    
    return translated_sequence


def correct_first_codon(data):
    base_tranlsate = nucleotide_to_protein(data=data, codons=3)
    

    print(base_tranlsate)
    
    mask_df_reduced = (data['protein_sequence'].str[0:7] != base_tranlsate)
    
    base_difference = data[mask_df_reduced]
    
    print(base_difference)


    '''
    if data['protein_sequence'][0] != first_base_tranlsate[0]:
        data['nucleotide_sequence'].str[0] = "A"
        '''
    
        
def main():
    # df_reduced = df_reduced.apply(correct_first_codon, axis=1)
    
    # print(df_reduced)

    correct_first_codon(df_reduced)

if __name__ == '__main__':
    main()