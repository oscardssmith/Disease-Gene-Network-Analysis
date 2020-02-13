import csv

LOOKUP_TABLE_PATH = "Data/protein-name-lookup-table.tsv"

def load_lookup_table():
    table = {}
    with open(LOOKUP_TABLE_PATH, newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)
        for row in reader:
            table[row[2]] = row[1]
    return table


def name_to_string(table, n):
    try:
        return list(table.keys())[list(table.values()).index(n)]
    except KeyError:
        return n


def string_to_name(table, s):
    try:
        return table[s]
    except KeyError:
        return s
    
