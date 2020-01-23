import networkx as nx
import numpy as np


def load_graph(path):
    """
    Loads data from TSV file pointed to by path into a networkx graph
    """
    graph = nx.Graph()
    with open(path, 'r') as input_file:
        input_file.readline()
        for line in input_file:
            # if line.strip() == "#":
            #    break
            data = line.strip().split(" ")
            #print("new edge:", data[0], data[1], data[2])
            graph.add_edge(data[0], data[1], confidence=data[2])
    return graph


def load_test_graph():
    return nx.complete_graph(5)


def load_test_start_vector():
    return np.asarray([0, 0, 0, 1, 0])


def load_disease_genes(path):
    """
    Loads disease genes from TSV file and returns a python list of all names
    """
    with open(path, 'r') as input_file:
        return input_file.read().splitlines()


def load_start_vector(path, ppi_graph):
    disease_gene_list = set(load_disease_genes(path))
    start_vector = np.zeros(ppi_graph.number_of_nodes())
    num_disease_genes = len(disease_gene_list)
    for i, node in enumerate(ppi_graph.nodes()):
        if node in disease_gene_list:
            start_vector[i] = 1/num_disease_genes
    return np.asarray(start_vector, order='F')
