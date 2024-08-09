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
from requests.adapters import HTTPAdapter, Retry


class InterProFetcher:
    def __init__(self):
        self.count = 0
        
    async def fetch(self, url, session):
        """Fetches a url, using specified ClientSession."""
        # using aysyncronous session to fetch data from the Interpro API
        async with session.get(url) as response:
            # check if the response is successful
            response.raise_for_status() 
            
            # add a count to keep track of the number of pages fetched
            self.count += 1
            print(f"Fetching page {self.count}...", end="\r")
            
            # await the response and return the json data
            return await response.json()

    async def fetch_all(self, url, session):
        """Fetches all pages of results"""
        # fetch all pages of results from the Interpro API
        results = []
        while url:
            # await the fetch function and extend the results list
            data = await self.fetch(url, session)
            results.extend(data["results"])
            # check if there is a next page
            url = data.get("next") 
        return results
    
    async def run(self, family_id, output_path):
        """Fetch the family members from the InterPro API asynchronously"""
        base_url = f"https://www.ebi.ac.uk/interpro/api/protein/uniprot/entry/interpro/{family_id}?page_size=200"
        
        async with ClientSession() as session:
            results = await self.fetch_all(base_url, session)
                
            with open(f"{output_path}/{family_id}_fetched.txt", "w") as file:
                accessions = [data["metadata"]["accession"] for data in results]
                file.write(",".join(accessions))
            print(f"Accessions have been written to {output_path}")

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
    
    # create an instance of the InterProFetcher class
    fetcher = InterProFetcher()
    
    # run the fetcher asynchronously
    await fetcher.run(args.family_id, args.output_path)
            
if __name__ == "__main__": 
    asyncio.run(main())