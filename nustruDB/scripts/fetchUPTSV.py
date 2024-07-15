#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script is adapted from https://www.uniprot.org/help/api_queries


usage: fetchUPTSV.py [-h] -i INPUT_FILE -o OUTPUT_PATH -n NAME [-w]

Retrieve uniprot features from uniprot IDs and store them in a tsv file

options:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input INPUT_FILE
                        Input file with uniprot IDs.
  -o OUTPUT_PATH, --output OUTPUT_PATH
                        Output path or directory to store the data.
  -n NAME, --name NAME  Name of the output file.
  -w, --overwrite       If file name already exists, overwrite it. Default is False.
  
Uniprots ids from a file or from the interpro script are used to fetch uniprot features from the Uniprot API and store them in a tsv file.
Features include: accession, gene primary, organism name, subcellular location, sequence, embl full, ccds, refseq full, alternative products, and sequence variation.
Example: python fetchUPTSV.py -i Example/examples_nustruDB/example_uniprotIDs.txt -o uniprot_features -n uniprot_features -w
"""

import re
import time
import json
import zlib

import argparse
import pandas as pd

from pathlib import Path

import requests
from urllib.parse import urlparse, parse_qs, urlencode
from requests.adapters import HTTPAdapter, Retry

# set the polling interval to 3 seconds and the API URL to the Uniprot API
POLLING_INTERVAL = 3
API_URL = "https://rest.uniprot.org"

# create a session with retries, with a total of 5 retries, a backoff factor of 0.25, and status codes 500, 502, 503, and 504
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries)) # mount the session with the retries

def check_response(response):
    """Checks if the response is successful, otherwise raise an exception."""
    try:
        response.raise_for_status() # debug requests response for errors
    # if the response is not successful, print the response and raise an exception
    except requests.HTTPError:
        print(response.json())
        raise

def submit_id_mapping(from_db, to_db, ids):
    """Submits a job to the ID mapping of Uniprot and return the job ID"""
    request = requests.post(
        f"{API_URL}/idmapping/run", # post the request to the Uniprot API and run the idmapping
        data={"from": from_db, "to": to_db, "ids": ",".join(ids)}, # data to be mapped
    )
    # check if the response is successful
    check_response(request)
    
    return request.json()["jobId"] # return the job ID

def get_next_link(headers):
    """Gets the next link from the headers of the response."""
    re_next_link = re.compile(r'<(.+)>; rel="next"') # use regex to find the next link
    # check if the link is in the headers and return the next link
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def check_id_mapping_results_ready(job_id):
    """Checks if the ID mapping results are ready"""
    while True: # loop until the results are ready
        # get the status of the job ID and check if the job is running
        request = session.get(f"{API_URL}/idmapping/status/{job_id}")
        check_response(request)
        j = request.json()
        # if the job is running, print a message and retry after the polling interval
        if "jobStatus" in j:
            if j["jobStatus"] == "RUNNING":
                print(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                # if the job is not running, raise an exception
                raise Exception(j["jobStatus"])
        else:
            # if the job is not running, return the results
            return bool(j["results"] or j["failedIds"])

def get_batch(batch_response, file_format, compressed):
    """Gets the batch results from the response"""
    # get the next link from the response of the batch
    batch_url = get_next_link(batch_response.headers)
    while batch_url: # loop until no next batch URL is available
        batch_response = session.get(batch_url)
        batch_response.raise_for_status()
        # decode the results from the batch response
        # return the decoded results
        yield decode_results(batch_response, file_format, compressed)
        batch_url = get_next_link(batch_response.headers)

def combine_batches(all_results, batch_results, file_format):
    """Combines the batch results with all results"""
    # add the results from the batch to the all results depending on the file format
    if file_format == "json": 
        # check keys in the batch results and add them to the all results
        for key in ("results", "failedIds"):
            if key in batch_results and batch_results[key]:
                all_results[key] += batch_results[key]
    # return the all results with the batch results if available
    elif file_format == "tsv":
        return all_results + batch_results[1:]
    else:
        return all_results + batch_results
    return all_results

def get_id_mapping_results_link(job_id):
    """Gets the ID mapping results link"""
    # get the redirect URL of the job ID and return the redirect URL
    url = f"{API_URL}/idmapping/details/{job_id}"
    request = session.get(url)
    check_response(request)
    return request.json()["redirectURL"]

def decode_results(response, file_format, compressed):
    """Decodes the results from the response"""
    # decode (decompress) the results from the response depending on the file format
    if compressed:
        # use zlib to decompress the response content
        decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
        # depending on the file format, decode the decompressed content wih decoding
        if file_format == "json":
            j = json.loads(decompressed.decode("utf-8"))
            return j
        elif file_format == "tsv":
            # return the decompressed content as a list of lines by splitting the content by new line
            return [line for line in decompressed.decode("utf-8").split("\n") if line]
        elif file_format == "xlsx":
            return [decompressed]
        elif file_format == "xml":
            return [decompressed.decode("utf-8")]
        else:
            return decompressed.decode("utf-8")
    # if the response is not compressed, return the response content depending on the file format
    elif file_format == "json":
        return response.json()
    elif file_format == "tsv":
        return [line for line in response.text.split("\n") if line]
    elif file_format == "xlsx":
        return [response.content]
    elif file_format == "xml":
        return [response.text]
    return response.text

def print_progress_batches(batch_index, size, total):
    """Prints the progress of the batches"""    
    # calculate the number of fetched results and print the progress
    n_fetched = min((batch_index + 1) * size, total)
    print(f"Fetched: {n_fetched} / {total}")

def get_id_mapping_results_search(url):
    """Gets the ID mapping results from the search URL"""
    parsed = urlparse(url) # parse the URL
    query = parse_qs(parsed.query) # parse the query string
    
    # get the file format, size, and compressed from the query string
    # set the default file format to json and size to 500 if not defined in the query
    # also check if the results should be compressed or not
    file_format = query["format"][0] if "format" in query else "json"
    if "size" in query:
        size = int(query["size"][0])
    else:
        size = 500
        query["size"] = size
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    
    # replace the query string with the updated query string
    parsed = parsed._replace(query=urlencode(query, doseq=True))
    url = parsed.geturl() # get the updated URL
    request = session.get(url) # get the request from the updated URL
    check_response(request)
    
    # assign results to the decoded results from the request in the specified file format
    # assign total to the total results from the request headers
    results = decode_results(request, file_format, compressed)
    total = int(request.headers["x-total-results"])
    
    
    print_progress_batches(0, size, total) # print the progress of the batches
    # loop through the batches and combine the results
    for i, batch in enumerate(get_batch(request, file_format, compressed), 1):
        results = combine_batches(results, batch, file_format)
        print_progress_batches(i, size, total) # update the progress of the batches fetched
    return results

def main() -> None:
    parser = argparse.ArgumentParser(
        prog='fetchUPTSV.py',
        description="Retrieve uniprot features from uniprot IDs and store them in a tsv file"
    )
    parser.add_argument(
        '-i', '--input', type=str, dest="input_file", required=True,
        help='Input file with uniprot IDs.'
    )
    parser.add_argument( 
        '-o', '--output', type=str, dest="output_path", required=True,
        help='Output path or directory to store the data.'
    )
    parser.add_argument(
        '-n', '--name', type=str, dest="name", required=True,
        help='Name of the output file.'
    )
    parser.add_argument(
        '-w', '--overwrite', action="store_true", dest="overwrite", required=False, default=False,
        help='If file name already exists, overwrite it. Default is False.' 
    )
    args = parser.parse_args()
    
    # read the uniprot IDs from the input file
    with open(args.input_file, "r") as f:
        ids = f.read().splitlines()
        
    # submit the ID mapping job and get the job ID
    job_id = submit_id_mapping(from_db="UniProtKB_AC-ID", to_db="UniProtKB", ids=ids)
    
    # check if the ID mapping results are ready and get the results
    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        # get the ID mapping results from the search URL and store them in a tsv file with the defined features as columns
        results = get_id_mapping_results_search(link + "?format=tsv" + "&fields=accession%2Cgene_primary%2Corganism_name%2Ccc_subcellular_location%2Csequence%2Cxref_embl_full%2Cxref_ccds%2Cxref_refseq_full%2Ccc_alternative_products%2Cft_var_seq") 

    # store the results in a tsv file with the defined features as columns 
    columns = results[0].split("\t") 
    data = [line.split("\t") for line in results[1:]]
    
    df = pd.DataFrame(data, columns=columns) # create a dataframe from the data
    df.drop(columns=["From"], inplace=True) # remove the "From" column
    
    # check if the output path exists and write the dataframe to a tsv file
    if Path(args.output_path).exists() or args.overwrite:
        df.to_csv(f"{args.output_path}/{args.name}.tsv", sep="\t", index=False)
    else:
        print(f"Error: {args.output_path} already exists. Use -w to overwrite.")
        exit(1)
    
if __name__ == "__main__":
    main()