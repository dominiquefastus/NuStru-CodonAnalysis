#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio import AlignIO
from Bio import SeqIO
import python_codon_tables as pct
import pandas as pd
import numpy as np


from plotly.subplots import make_subplots
import plotly.graph_objects as go


name_of_protein_alignment = "rtx_prot_fam_ecoli_protein_aligned.fasta"
name_of_nucleotide = "rtx_prot_fam_ecoli_nucleotide.fasta"

working_dir = "/Users/dominiquefastus/Downloads/nustruTREE/MSA"

protein_alignment = AlignIO.read(f"{working_dir}/{name_of_protein_alignment}", "fasta")
nucleotide_alignment = SeqIO.parse(f"{working_dir}/{name_of_nucleotide}", "fasta")

nustrudb = pd.read_csv("/Users/dominiquefastus/Downloads/rtx_prot_fam_ecoli_secstru.csv")

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

def cub_msa_table(prot_seq_arr=None, cod_seq_arr=None):
    cub_table = {
    # '*': {'TAA': None, 'TAG': None, 'TGA': None}, ignoring stop codons
    'A': {'GCA': None, 'GCC': None, 'GCG': None, 'GCT': None},
    'C': {'TGC': None, 'TGT': None},
    'D': {'GAC': None, 'GAT': None},
    'E': {'GAA': None, 'GAG': None},
    'F': {'TTC': None, 'TTT': None},
    'G': {'GGA': None, 'GGC': None, 'GGG': None, 'GGT': None},
    'H': {'CAC': None, 'CAT': None},
    'I': {'ATA': None, 'ATC': None, 'ATT': None},
    'K': {'AAA': None, 'AAG': None},
    'L': {'CTA': None, 'CTC': None, 'CTG': None, 'CTT': None, 'TTA': None, 'TTG': None},
    'M': {'ATG': None},
    'N': {'AAC': None, 'AAT': None},
    'P': {'CCA': None, 'CCC': None, 'CCG': None, 'CCT': None},
    'Q': {'CAA': None, 'CAG': None},
    'R': {'AGA': None, 'AGG': None, 'CGA': None, 'CGC': None, 'CGG': None, 'CGT': None},
    'S': {'AGC': None, 'AGT': None, 'TCA': None, 'TCC': None, 'TCG': None, 'TCT': None},
    'T': {'ACA': None, 'ACC': None, 'ACG': None, 'ACT': None},
    'V': {'GTA': None, 'GTC': None, 'GTG': None, 'GTT': None},
    'W': {'TGG': None},
    'Y': {'TAC': None, 'TAT': None}}
    
    for aa in cub_table.keys():
        n_AA = np.count_nonzero(prot_seq_arr == aa)
                
        nc_AA = len(cub_table[aa].keys())
        
        for codon in cub_table[aa].keys():
            nc = np.count_nonzero(cod_seq_arr == codon)
            
            fc =(nc / n_AA) * 1/nc_AA

            cub_table[aa][codon] = round(fc,6)
        
    return cub_table

def map_rarity(protein_alignment, nustrudb, cu_table):
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
                alignment_value_matrix[i, position] = cu_table[aa][sequence[position_adj*3:position_adj*3+3].upper()]

    residue_mean = []
    for col_mean in np.sum(alignment_value_matrix, axis=0):
        residue_mean.append(col_mean / len(seq_name))
    
    return alignment_value_matrix, seq_name, seq_pos, residue_mean

all_seqs_protein = fasta_to_array(protein_alignment, codon=False)
all_seqs_nt = fasta_to_array(fasta=nucleotide_alignment, align_to=all_seqs_protein, codon=True)

cub_msa_table_ddla = cub_msa_table(prot_seq_arr=all_seqs_protein, cod_seq_arr=all_seqs_nt)
alignment_value_matrix, seq_name, seq_pos, residue_mean = map_rarity(protein_alignment, nustrudb, cub_msa_table_ddla)
sorted_alignment_value_matrix = np.sort(alignment_value_matrix, axis=0)
sorted_alignment_value_matrix = np.flip(sorted_alignment_value_matrix, axis=0)

fig = make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.02)
fig.add_trace((go.Heatmap(z=sorted_alignment_value_matrix, y=seq_name, colorscale="reds")), row=3, col=1)
fig.add_trace((go.Heatmap(z=alignment_value_matrix, y=seq_name, colorscale="blues")), row=2, col=1)
fig.add_trace((go.Scatter(x=seq_pos, y=residue_mean, mode="lines", fill="toself", line=dict(color="royalblue"))), row=1, col=1)
fig.show()
