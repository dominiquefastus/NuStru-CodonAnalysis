#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio import AlignIO
from Bio import SeqIO
import python_codon_tables as pct
import pandas as pd
import numpy as np
import pathlib
import math
import os 

from plotly.subplots import make_subplots
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import seaborn as sns

name_of_protein_alignment = "test_ddla_protein_aligned.fasta"
name_of_nucleotide_alignment = "test_ddla_nucleotide_aligned.fasta"

working_dir = "/Users/dominiquefastus/master_project/NuStru/nustruEVOL/nustruTREE/MSA"

protein_alignment = AlignIO.read(f"{working_dir}/{name_of_protein_alignment}", "fasta")
nucleotide_alignment = SeqIO.parse(f"{working_dir}/{name_of_nucleotide_alignment}", "fasta")
nustrudb = pd.read_csv("/Users/dominiquefastus/master_project/NuStru/nustruDB/DDLA_uniprot_sec_struct_04.csv")

e_coli_pct = pct.get_codons_table("e_coli_316407")

def fasta_to_array(fasta, align_to=None, codon=False):
    all_seqs = []
    all_ids = []
    
    if codon:
        for (ind,record) in enumerate(fasta):
            all_seqs.append(list(str(record.seq)))
            all_seqs[ind] = [''.join(map(str, all_seqs[ind][i:i+3])) for i in range(0, len(all_seqs[ind]), 3)]
            
    else:           
        for record in fasta:
            all_seqs.append(list(str(record.seq)))
            all_ids.append(record.id)
        
        all_seqs = np.array(all_seqs)

    if align_to is not None:
        gap_indeces = np.where(align_to == '-')
        
        for gap_index in zip(gap_indeces[0], gap_indeces[1]):
            all_seqs[gap_index[0]].insert(gap_index[1], '---')
        
        all_seqs = np.array(all_seqs)
        # deleting stop codons as no protein assigned to them
        all_seqs = np.delete(all_seqs, -1, axis=1)
            
    
    
    all_ids = np.array(all_ids).reshape(len(all_ids), 1)
    # all_seqs = np.append(all_seqs, all_ids, axis=1)
    
    return all_seqs


def map_rarity(protein_alignment, nustrudb, e_coli_pct):
    codon_position_start = 0
    alignment_value_matrix = np.zeros((len(protein_alignment), len(protein_alignment[0])))
    seq_name = [seq.id for seq in protein_alignment]
    seq_pos = [i for i in range(len(protein_alignment[0]))]

    sart_count = [0 for i in range(len(seq_name))]
    pos_count_dict = {seq_name[i]: 0 for i in range(len(seq_name))}
    for position in range(len(protein_alignment[0])):

        for i, (aa, seq) in enumerate(zip(protein_alignment[:,position],seq_name)):
            if aa == '-':
                alignment_value_matrix[i, position] = 0
                pos_count_dict[seq] += 1
            else:
                prot_position = pos_count_dict[seq]
                position_adj = position - prot_position
                
                sequence = nustrudb[nustrudb["primary_id"] == seq]["nucleotide_sequence"].values[0]
                alignment_value_matrix[i, position] = e_coli_pct[aa][sequence[position_adj*3:position_adj*3+3].upper()]

    residue_mean = []
    for col_mean in np.mean(alignment_value_matrix, axis=0):
        residue_mean.append(col_mean)
        
    return alignment_value_matrix, seq_name, seq_pos, residue_mean

alignment_value_matrix, seq_name, seq_pos, residue_mean = map_rarity(protein_alignment, nustrudb, e_coli_pct)


fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.02)
fig.add_trace((go.Heatmap(z=alignment_value_matrix, y=seq_name, colorscale="blues")), row=2, col=1)
fig.add_trace((go.Bar(x=seq_pos, y=residue_mean)), row=1, col=1)
fig.show()