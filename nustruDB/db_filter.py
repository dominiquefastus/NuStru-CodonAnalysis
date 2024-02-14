#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq

df = pd.read_csv('/Users/dominiquefastus/master_project/NuStru/nustruDB/E_COLI_K12_02.csv', index_col=False)


df_reduced = df.drop_duplicates(subset=['primary_id','nucleotide_sequence'])

def nucleotide_to_protein(data, codons = 3):
    translated_sequence = data['nucleotide_sequence'].str[:codons].apply(lambda x: Seq(x).translate())
    
    return translated_sequence


def check_first_codon(data):
    first_base_tranlsate = nucleotide_to_protein(data=data, codons=3)

    if data['protein_sequence'][0] != first_base_tranlsate:
        print(data['protein_sequence'][0])

def main():
    check_first_codon(df_reduced)

if __name__ == '__main__':
    print('test')