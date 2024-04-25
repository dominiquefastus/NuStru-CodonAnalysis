#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio import AlignIO
from Bio import SeqIO
import pandas as pd
import numpy as np
import argparse
import os

from ete3 import Tree

from plotly.subplots import make_subplots
import plotly.graph_objects as go

import matplotlib.pyplot as plt
import seaborn as sns

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

    residue_sum = []
    for col_mean in np.sum(alignment_value_matrix, axis=0):
        residue_sum.append(col_mean / len(seq_name))
        
    residue_max = []
    for col_max in np.max(alignment_value_matrix, axis=0):
        residue_max.append(col_max)
        
    sequence_sum = {}
    for col_mean, name in zip(np.sum(alignment_value_matrix, axis=1), seq_name):
        sequence_sum[name] = round((col_mean / len(alignment_value_matrix[0])), 5)
    
    return alignment_value_matrix, seq_name, seq_pos, residue_sum, residue_max, sequence_sum

def sum_distances_to_root(tree):
    leaf_distances = {}
    for leaf in tree.iter_leaves():
        current, distance_sum = leaf, 0
        while not current.is_root():
            distance_sum += current.dist
            current = current.up
        leaf_distances[leaf.name] = round(distance_sum, 5)
    return leaf_distances

def main() -> None:
    
    parser = argparse.ArgumentParser(
        prog='msa_codon_rarity.py',
        description="Calculate the codon rarity for a given protein alignment and visualize the results."
    )
    parser.add_argument(
        '-a', '--alignment', type=str, dest="alignment_file", required=True,
        help='Input file of aligned protein sequences in fasta format.'
    )
    parser.add_argument(
        '-nt', '--nucleotide', type=str, dest="nt_fasta_file", required=True,
        help='''Input file of corresponding nucleotide sequences of the aligned protein sequence in fasta format. 
                The nucleotide don't have to be aligned or ordered, but shoud have the same sequence identifiers as the protein alignment.'''
    )
    parser.add_argument(
        '-db', '--database', type=str, dest="nustruDB", required=True,
        help='''Database (like nustruDB) to map the protein sequence to the nucleotide sequence and secondary structure. 
                Database should be in csv format and include the primary_id, the protein and nucleotide sequence, and secondary structure.'''
    )
    parser.add_argument(
        '-t', '--tree', type=str, dest="tree_file", required=False, default=None,
        help='Constructed tree of alignment in newick format to calculate the divergence of the codon rarity.'
    )
    parser.add_argument( 
        '-o', '--output', type=str, dest="output_path", default=os.getcwd(),
        help='Output to store the plots and results.'
    )
    
    args = parser.parse_args()
    
    protein_alignment = AlignIO.read(args.alignment_file, "fasta")
    nucleotide_sequences = SeqIO.parse(args.nt_fasta_file, "fasta")
    nustrudb = pd.read_csv(args.nustruDB)
    
    all_seqs_protein = fasta_to_array(protein_alignment, codon=False)
    all_seqs_nt = fasta_to_array(fasta=nucleotide_sequences, align_to=all_seqs_protein, codon=True)
    
    cub_msa_table_ddla = cub_msa_table(prot_seq_arr=all_seqs_protein, cod_seq_arr=all_seqs_nt)
    alignment_value_matrix, seq_name, seq_pos, residue_sum, residue_max, sequence_sum = map_rarity(protein_alignment, nustrudb, cub_msa_table_ddla)
    
    sorted_alignment_value_matrix = np.sort(alignment_value_matrix, axis=0)
    sorted_alignment_value_matrix = np.flip(sorted_alignment_value_matrix, axis=0)

    fig = make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.02)
    fig.add_trace((go.Heatmap(z=sorted_alignment_value_matrix, y=seq_name, colorscale="reds")), row=3, col=1)
    fig.add_trace((go.Heatmap(z=alignment_value_matrix, y=seq_name, colorscale="blues")), row=2, col=1)
    fig.add_trace((go.Scatter(x=seq_pos, y=residue_max, mode="lines", fill="toself", line=dict(color="lightcoral"))), row=1, col=1)
    fig.write_image(f"{args.output_path}/codon_rarity_heatmap.png")
    fig.write_html(f"{args.output_path}/codon_rarity_heatmap.html")
    
    
    if args.tree_file is not None:
        tree = Tree(args.tree_file)
        leaf_distances = sum_distances_to_root(tree)
        x_values = [leaf_distances[key] for key in sorted(leaf_distances)]
        y_values = [sequence_sum[key] for key in sorted(sequence_sum)]
        
        '''
        keys = sorted(sequence_sum.keys())
        
        # Optionally, label each point
        for i, key in enumerate(keys):
            plt.text(x_values[i], y_values[i], key)
        '''
        
        plt.figure(figsize=(8, 5))
        plt.style.use('ggplot')
        sns.regplot(x=x_values, y=y_values,
                    scatter_kws={'color': 'blue', 'alpha': 0.5}, 
                    line_kws={'color': 'navy', })
        plt.title('Divergence of Codon Rarity')
        plt.xlabel('Branch Lengths from Root (divergence)')
        plt.ylabel('Codon Rarity Score (sum of all residues)')
        plt.savefig(f"{args.output_path}/divergence_codon_rarity.png")
        

if __name__ == "__main__":
    main()