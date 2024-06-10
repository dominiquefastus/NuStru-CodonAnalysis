#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: fetchINPRO.py [-h] -o OUTPUT_PATH [-w] family_id

Fetch Accessions from InterPro

positional arguments:
  family_id             InterPro family ID to fetch accessions for

options:
  -h, --help            show this help message and exit
  -o OUTPUT_PATH, --output OUTPUT_PATH
                        Output to store the interpro accessions.
  -w, --overwrite       If file name already exists, overwrite it. Default is False.
  
Use an InterPro family ID to fetch Uniprot and PDB accessions (protein family members from the InterPro API and store them in a file.
Example: python fetchINPRO.py IPR000839 -o IPR000839_accessions.txt
"""

import asyncio
import argparse
from pathlib import Path
from aiohttp import ClientSession

count = 0
async def fetch(url, session):
    """Fetches a url, using specified ClientSession"""
    # using aysyncronous session to fetch data from the Interpro API
    async with session.get(url) as response:
        # check if the response is successful
        response.raise_for_status() 
        
        # add a count to keep track of the number of pages fetched
        global count
        count += 1
        print(f"Fetching page {count}...", end="\r")
        
        # await the response and return the json data
        return await response.json()

async def fetch_all(url, session):
    """Fetches all pages of results"""
    # fetch all pages of results from the Interpro API
    results = []
    while url:
        # await the fetch function and extend the results list
        data = await fetch(url, session)
        results.extend(data["results"])
        # check if there is a next page
        url = data.get("next") 
    return results

async def main():
    parser = argparse.ArgumentParser(
        prog="fetchINPRO.py",
        description="Fetch Accessions from InterPro"
    )
    parser.add_argument(
        "family_id", type=str,
        help="InterPro family ID to fetch accessions for"
    )
    parser.add_argument( 
        '-o', '--output', type=str, dest="output_path", required=True,
        help='Output to store the interpro accessions.'
    )
    parser.add_argument(
        '-w', '--overwrite', action="store_true", dest="overwrite", required=False, default=False,
        help='If file name already exists, overwrite it. Default is False.' 
    )
    args = parser.parse_args()
    
    # set the base url to fetch the data from the Interpro API (200 results per page)
    base_url = f"https://www.ebi.ac.uk/interpro/api/protein/uniprot/entry/interpro/{args.family_id}?page_size=200"

    # using the aiohttp library to create a ClientSession
    async with ClientSession() as session:
        # await the fetch_all function to get all the results
        results = await fetch_all(base_url, session)
        if not Path(args.output_path).exists() or args.overwrite:
            with open(args.output_path, "w") as file: # write the accessions to a file
                for item in results:
                    accession = item["metadata"]["accession"]
                    file.write(f"{accession},")
            print(f"Accessions have been written to {args.output_path}")
        else:
            print(f"Error: {args.output_path} already exists. Use -w to overwrite.")
            exit(1)
            
if __name__ == "__main__": 
    asyncio.run(main())