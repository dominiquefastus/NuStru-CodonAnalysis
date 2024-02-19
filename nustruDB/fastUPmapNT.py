import pandas as pd 

from Bio import Entrez
from typing import List

def insert_pandas(df, source, entry_id, gene_name, organism, expression_system, mitochondrial, protein_sequence, nucleotide_id, nucleotide_sequence, plddt):
    df = df._append({"source": source, "primary_id": entry_id, "gene_name": gene_name, "organism": organism, "expression_system": expression_system, "mitochondrial": mitochondrial,
                     "protein_sequence": protein_sequence, "nucleotide_id": nucleotide_id, "nucleotide_sequence": nucleotide_sequence, "plddt": plddt}, ignore_index=True)
    
    # print(f'Entry {entry_id} successfully inserted in pandas DataFrame!')
    
    return df

def retrieve_nucleotide_seq(ncbiDB="nuccore", entryID=None, rettype="fasta", retmode="text"):
    Entrez.email = "dominique.fastus@biochemistry.lu.se"    
    handle = Entrez.efetch(db=ncbiDB, id=entryID, rettype=rettype, retmode=retmode)
    
    return handle.read()

def filter_sequence(sequences: str, searchID: List) -> None:
    search_criteria = searchID
    unique_sequences = set()
    matching_data = []
    entries = sequences.split('>')

    for entry in entries:
        if any(criterion in entry for criterion in search_criteria):
            header, sequence = entry.split('\n', 1)
            sequence = sequence.strip()
            # Check if the sequence is already added to ensure uniqueness
            if sequence not in unique_sequences:
                unique_sequences.add(sequence)
                matching_data.append((header, sequence.replace('\n',"")))

    return matching_data

def get_cds(data):
    id = data['nucleotide_id'].strip().replace('"', '').replace(' ', '')
    id_list = id.split(';')
    del id_list[-1]
    
    cds_ids = id_list[0::4]
    protein_ids = id_list[1::4]
    status = id_list[2::4]
    
    for cds_id, protein_id, status in zip(cds_ids, protein_ids, status):
        if status == '-':
            if retrieve_nucleotide_seq(ncbiDB="protein", entryID=protein_id).partition("\n")[2].strip().replace('\n','') == data['protein_sequence']:
                response_sequences = retrieve_nucleotide_seq(entryID=cds_id, rettype="fasta_cds_na")
                
                matched_id, matched_seq = filter_sequence(sequences=response_sequences,searchID=[data['primary_id'], protein_id])[0]  

                    
                break  
                    
    with open('/Users/dominiquefastus/Downloads/sequences.fasta', 'a') as f:
        f.write(f'>{matched_id}\n{matched_seq}\n')
    
    return matched_seq

df = pd.read_csv('/Users/dominiquefastus/Downloads/uniprotkb_taxonomy_id_562_2024_02_19.tsv', sep='\t', header=0, nrows=2)


nucleotide_protein_seqs_df = df[['Entry', 'Gene Names (primary)', 'Organism', 'Subcellular location [CC]', 'Sequence', 'EMBL']]
nucleotide_protein_seqs_df.rename(columns={'Entry': 'primary_id', 'Gene Names (primary)': 'gene_name', 'Organism': 'organism', 'Subcellular location [CC]': 'expression_system', 'Sequence': 'protein_sequence', 'EMBL': 'nucleotide_id'}, inplace=True)

print(nucleotide_protein_seqs_df.apply(get_cds, axis=1))