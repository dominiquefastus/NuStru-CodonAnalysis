import aiohttp
import asyncio
import csv
import argparse
from io import StringIO

async def fetch_data(session, url):
    async with session.get(url) as response:
        return await response.text()

async def download_uniprot_data(url, query_params, output_file):
    async with aiohttp.ClientSession() as session:
        page = 0
        with open(output_file, 'w', newline='') as file:
            writer = None
            while True:
                query_params['offset'] = page * query_params['size']
                full_url = f"{url}?{'&'.join(f'{k}={v}' for k, v in query_params.items())}"
                data = await fetch_data(session, full_url)

                reader = csv.DictReader(StringIO(data), delimiter='\t')
                records = list(reader)
                if not records:
                    break  

                if not writer:
                    writer = csv.DictWriter(file, fieldnames=reader.fieldnames, delimiter='\t')
                    writer.writeheader()
                writer.writerows(records)

                page += 1
                if len(records) < query_params['size']:
                    break  
                
async def main():
    parser = argparse.ArgumentParser(description="Download data from UniProt by given ID and write to TSV file.")
    parser.add_argument("uniprotid", type=str, help="UniProt ID to query")
    parser.add_argument("--output", type=str, help="Path to output TSV file")
    args = parser.parse_args()
    
    url = 'https://rest.uniprot.org/uniprotkb/search'
    query_params = {
        'fields': 'accession,gene_primary,organism_name,cc_subcellular_location,sequence,xref_embl_full,xref_ccds,xref_refseq_full,cc_alternative_products,ft_var_seq',
        'format': 'tsv',
        'query': f'({args.uniprotid})',
        'size': 500
    }
    await download_uniprot_data(url, query_params, args.output)

if __name__ == "__main__":
    asyncio.run(main())
