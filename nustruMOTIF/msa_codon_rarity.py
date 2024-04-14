
from Bio import AlignIO
import matplotlib.pyplot as plt
import python_codon_tables as pct
import seaborn as sns
import pandas as pd
import numpy as np
import pathlib
import os 

name_of_protein_alignment = "test_ddla_protein_aligned.fasta"
name_of_nucleotide_alignment = "test_ddla_nucleotide_aligned.fasta"

working_dir = "/Users/dominiquefastus/master_project/NuStru/nustruEVOL/nustruTREE/MSA"

protein_alignment = AlignIO.read(f"{working_dir}/{name_of_protein_alignment}", "fasta")
nucleotide_alignment = AlignIO.read(f"{working_dir}/{name_of_nucleotide_alignment}", "fasta")
nustrudb = pd.read_csv("/Users/dominiquefastus/master_project/NuStru/nustruDB/DDLA_uniprot_sec_struct_04.csv")

e_coli_pct = pct.get_codons_table("e_coli_316407")

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

import plotly.express as px
import plotly.graph_objects as go

from plotly.subplots import make_subplots
import plotly.graph_objects as go


fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.02)
fig.add_trace((go.Heatmap(z=alignment_value_matrix, y=seq_name)), row=2, col=1)
fig.add_trace((go.Bar(x=seq_pos, y=residue_mean)), row=1, col=1)
fig.show()