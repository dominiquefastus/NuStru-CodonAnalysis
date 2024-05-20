#!/usr/bin/env python
# -*- coding: utf-8 -*-

import asyncio
import argparse
from pathlib import Path
from aiohttp import ClientSession

count = 0
async def fetch(url, session):
    """Fetch a url, using specified ClientSession."""
    async with session.get(url) as response:
        response.raise_for_status() 
        global count
        count += 1
        print(f"Fetching page {count}...", end="\r")
        return await response.json()

async def fetch_all(url, session):
    """Fetch all pages of results."""
    results = []
    while url:
        data = await fetch(url, session)
        results.extend(data["results"])
        url = data.get("next") 
    return results

async def main():
    parser = argparse.ArgumentParser(
        prog="fetchINPRO.py",
        description="Fetch Accessions from InterPro"
    )
    parser.add_argument(
        "family_id", required=True, type=str,
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
    
    base_url = f"https://www.ebi.ac.uk/interpro/api/protein/uniprot/entry/interpro/{args.family_id}?page_size=200"
    async with ClientSession() as session:
        results = await fetch_all(base_url, session)
        if Path(args.output_path).exists() or args.overwrite:
            with open(args.output_path, "w") as file:
                for item in results:
                    accession = item["metadata"]["accession"]
                    file.write(f"{accession},")
            print(f"Accessions have been written to {args.output}")
        else:
            print(f"Error: {args.output_path} already exists. Use -w to overwrite.")
            exit(1)
            

if __name__ == "__main__": 
    asyncio.run(main())
