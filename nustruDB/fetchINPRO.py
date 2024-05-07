import asyncio
import argparse
from aiohttp import ClientSession

async def fetch(url, session):
    """Fetch a url, using specified ClientSession."""
    async with session.get(url) as response:
        response.raise_for_status() 
        return await response.json()

async def fetch_all(url, session):
    """Fetch all pages of results."""
    results = []
    while url:
        data = await fetch(url, session)
        results.extend(data["results"])
        url = data.get("next") 
    return results

async def main(family_id, output_file):
    parser = argparse.ArgumentParser(description="Fetch Accessions from InterPro")
    parser.add_argument("family_id", help="InterPro family ID to fetch accessions for")
    parser.add_argument("--output", help="Output file to write accessions", default="accessions.txt")
    
    args = parser.parse_args()
    
    base_url = f"https://www.ebi.ac.uk/interpro/api/protein/uniprot/entry/interpro/{args.family_id}?page_size=200"
    async with ClientSession() as session:
        results = await fetch_all(base_url, session)
        with open(args.output_file, "w") as file:
            for item in results:
                accession = item["metadata"]["accession"]
                file.write(f"{accession}\n")
        print(f"Accessions have been written to {args.output_file}")

if __name__ == "__main__": 
    asyncio.run(main())
