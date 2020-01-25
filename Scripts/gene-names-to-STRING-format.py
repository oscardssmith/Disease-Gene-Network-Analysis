# Code based on (shamelessly taken and then modified) from
# https://string-db.org/cgi/help.pl?subpage=api%23mapping-identifiers
##########################################################
# For a given list of proteins the script resolves them
# (if possible) to the best matching STRING identifier
# and prints out the mapping on screen in the TSV format
###########################################################

import random
import requests
import sys

STRING_API_URL = "http://string-db.org/api"
OUTPUT_FORMAT = "tsv-no-header"
METHOD = "get_string_ids"

# File with common gene names
def read_genes(infile):
    genes = []
    with open(infile) as infile:
        for line in infile:
            # Don't read \n
            genes.append(line[:-1])
    print("num genes: ", len(genes))
    return genes

# contruct params dictionary
def gene_names_to_STRING(genes):
    params = {
        "identifiers": "\r".join(genes),  # your protein list
        "species": 9606,  # species NCBI identifier
        "limit": 1,  # only one (best) identifier per input protein
        "echo_query": 1,  # see your input identifiers in the output
        "caller_identity": "www.awesome_app.org"  # your app name
    }

    # contruct method URL
    request_url = STRING_API_URL + "/" + OUTPUT_FORMAT + "/" + METHOD

    # Call STRING
    try:
        response = requests.post(request_url, params=params)
    except requests.exceptions.RequestException as e:
        print(e)
        sys.exit()
    return response


def output_response(response):
    # Read and parse the results
    count_response = 0
    protein_file = open("protein_file.txt", "a")
    for line in response.text.strip().split("\n"):
        fields = line.split("\t")
        # print(fields)
        input_identifier_, string_identifier = fields[0], fields[2]
        #print("Input:", input_identifier, "STRING:", string_identifier, sep="\t")
        print(string_identifier)
        protein_file.write(string_identifier)
        protein_file.write("\n")
        count_response += 1
    print("num proteins outputted:", count_response)
    protein_file.close()


if __name__ == '__main__':
    if len(sys.argv) == 2:
        infile = sys.argv[1]
    genes = read_genes(infile)
    response = gene_names_to_STRING(genes)
    output_response(response)
