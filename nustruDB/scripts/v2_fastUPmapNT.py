#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
INFO: Pandarallel will run on 10 workers.
INFO: Pandarallel will use standard multiprocessing data transfer (pipe) to transfer data between the main process and workers.

usage: fastUPmapNT.py [-h] -i INPUT_FILE -o OUTPUT_PATH -n NAME [-w]

Retrieve nucleotide sequences from uniprot IDs in a csv file.

options:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input INPUT_FILE
                        Input file with uniprot IDs in a csv format.
  -o OUTPUT_PATH, --output OUTPUT_PATH
                        Output path or directory to store the log file and the data.
  -n NAME, --name NAME  Name of the output files and log file.
  -w, --overwrite       If file name already exists, overwrite it. Default is False.
  
Script to retrieve nucleotide sequences from uniprot IDs in a csv/ tsv file. The script will create a log file and a new csv file.
Example: python fastUPmapNT.py -i Example/examples_nustruDB/example_uniprotList.tsv -o . -n example_uniprotList_nustru [-w]
"""

import argparse
import aiohttp
import asyncio

import pandas as pd 
import numpy as np
pd.options.mode.chained_assignment = None
pd.set_option('future.no_silent_downcasting', True)

from Bio.Seq import Seq

from pathlib import Path
from itertools import chain
from typing import List
import logging
import sys

# set the number of workers for parallel processing to 10 (maximum for API requests limit per second)
# if the computer has less cores, all cores will be used
# from pandarallel import pandarallel
# pandarallel.initialize(progress_bar=True, nb_workers=10, verbose=0)

class UPmapperNT:
    def __init__(self, retmax=5_000, concurrency=100) -> None:
        self.retmax = retmax
        self.concurrency = concurrency
        self.fetched_nucleotide_seqs = []
        
        self.efetch_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                
    async def fetch_batch(self, session, current_ids, batch_num):
        logging.info(f"Batch {batch_num} fetching...")
        efetch_params = {
            "db": "protein",
            "id": ",".join(current_ids),
            "retmax": len(current_ids),
            "rettype": "fasta_cds_na",
            "retmode": "text"
        }
        async with session.post(self.efetch_base, data=efetch_params) as response:
            if response.status == 200:
                text = await response.text()
                logging.info(f"Batch {batch_num}: Fetched {len(current_ids)}.")
                seqs = text.replace("\n\n", "*").replace("\n","").replace("CDS]","*").split("*")[1::2]
                self.fetched_nucleotide_seqs.append(seqs)
            else:
                logging.error(f"Error in Batch {batch_num}: {response.status}")
                return None

    async def fetch_all(self, ids):
        async with aiohttp.ClientSession() as session:
            tasks = []
            batch_num = 1
            for i in range(0, len(ids), self.retmax):
                current_ids = ids[i:i+self.retmax]
                task = asyncio.create_task(self.fetch_batch(session, current_ids, batch_num))
                tasks.append(task)
                batch_num += 1
                await asyncio.sleep(1)
                if len(tasks) >= self.concurrency:
                    # Wait for current batch of tasks to complete before sending more requests
                    await asyncio.gather(*tasks)
                    tasks = []

            # If there are remaining tasks, ensure they finish
            if tasks:
                await asyncio.gather(*tasks)

    def merge_with_hits(self, nucleotide_protein_seqs_df, fetched_nucleotide_seqs) -> pd.DataFrame:
        def translate_nucleotide_seq(seq):
            # Translate only if length is divisible by 3
            try:
                if seq.endswith(("TAA","TAG","TGA")):
                    print("ok")
                    return Seq(seq).translate()[:-1]
                else:
                    remain = len(seq) % 3
                    return ''.join(Seq(seq[:-remain]).translate())
            except:
                return None
        
        # Apply translation to all fetched sequences
        translated_seqs = [translate_nucleotide_seq(seq) for seq in fetched_nucleotide_seqs]

        # Convert to dataframe for faster matching
        translated_df = pd.DataFrame({
            'nucleotide_sequence': fetched_nucleotide_seqs,
            'protein_sequence': translated_seqs
        })
        
        translated_df.to_csv("/Users/dominiquefastus/master_project/NuStru/Example/examples_nustruDB/tranlated_df.csv")

        nucleotide_protein_seqs_df = nucleotide_protein_seqs_df[['primary_id', 'gene_name', 'organism', 'expression_system', 'protein_sequence']]
        # Merge with original dataframe based on protein sequence match
        merged_df = nucleotide_protein_seqs_df.merge(translated_df, on='protein_sequence', how='inner')
        merged_df.dropna()
        
        return merged_df
    
    async def get_cds(self, nucleotide_protein_seqs_df, output_path, name):
        """Get CDS from uniprot IDs"""
        nucleotide_ids = [str(entry).strip().replace('"', '').replace(' ', '').split(";")[:-1] for entry in  nucleotide_protein_seqs_df["nucleotide_id"].to_list()]

        cds_ids = [entry[0::4] for entry in nucleotide_ids]
        cds_ids = list(chain.from_iterable(cds_ids))
        protein_ids = [entry[1::4] for entry in nucleotide_ids]
        protein_ids = list(chain.from_iterable(protein_ids))
        status = [entry[2::4] for entry in nucleotide_ids]
        status = list(chain.from_iterable(status))

        plain_cds_id = []
        plain_protein_id = []

        for cds_id, protein_id, status in zip(cds_ids, protein_ids, status):
            if status == "-" and protein_id != "-":
                plain_cds_id.append(cds_id)
                plain_protein_id.append(protein_id)
                
        ###### Run the ... ######
        await self.fetch_all(plain_protein_id)
        
        ###### Construct the finished dataframe with the nt sequences ######
        fetched_nucleotide_seqs = list(chain.from_iterable(self.fetched_nucleotide_seqs))
        fetched_nucleotide_seqs = list(set(fetched_nucleotide_seqs))
        merged_df = self.merge_with_hits(nucleotide_protein_seqs_df, fetched_nucleotide_seqs)
        
        merged_df.to_csv(f'{output_path}/{name}.csv', mode='a', index=False, header=False)
                    

async def main():
    parser = argparse.ArgumentParser(
        prog='fastUPmapNT.py',
        description="Retrieve nucleotide sequences from uniprot IDs in a tsv file."
    )
    parser.add_argument(
        '-i', '--input', type=str, dest="input_file", required=True,
        help='Input file with uniprot IDs in a csv format.'
    )
    parser.add_argument( 
        '-o', '--output', type=str, dest="output_path", required=True,
        help='Output path or directory to store the log file and the data.'
    )
    parser.add_argument(
        '-n', '--name', type=str, dest="name", required=True,
        help='Name of the output files and log file.'
    )
    parser.add_argument(
        '-w', '--overwrite', action="store_true", dest="overwrite", required=False, default=False,
        help='If file name already exists, overwrite it. Default is False.' 
    )

    args = parser.parse_args()
    
    ####### Create a log file ######
    logging.basicConfig(handlers=[logging.FileHandler(f'{args.output_path}/{args.name}.log'), logging.StreamHandler(sys.stdout)],
                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.DEBUG)
    
    ###### Read the input file and replace empty strings with NaN ######
    nucleotide_protein_seqs_df = pd.read_csv(args.input_file, sep='\t', header=0, dtype={"Subcellular location [CC]": object, "Alternative sequence": object})
    nucleotide_protein_seqs_df = nucleotide_protein_seqs_df.replace(r'^\s*$', np.nan, regex=True)

    nucleotide_protein_seqs_df = nucleotide_protein_seqs_df[['Entry', 'Gene Names (primary)', 'Organism', 'Subcellular location [CC]', 'Sequence', 'EMBL']]
    nucleotide_protein_seqs_df.columns = ['primary_id', 'gene_name', 'organism', 'expression_system', 'protein_sequence', 'nucleotide_id']
    
    ###### Create a csv file with the header (columns) ######
    if not Path(f'{args.output_path}/{args.name}.csv').exists() or args.overwrite:
        with open(f'{args.output_path}/{args.name}.csv', mode='w') as f:
            f.write('source,primary_id,gene_name,organism,expression_system,protein_sequence,nucleotide_sequence\n')
    else:
        logging.error(f"Error: {args.output_path}/{args.name}.csv already exists.")
        exit(1)
        
    ###### Create an instance of the UPmapperNT class and run ######
    mapper = UPmapperNT()
    await mapper.get_cds(nucleotide_protein_seqs_df=nucleotide_protein_seqs_df, 
                         output_path=args.output_path, name=args.name)
    
if __name__ == '__main__':
    asyncio.run(main())