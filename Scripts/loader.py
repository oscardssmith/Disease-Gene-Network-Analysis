import networkx as nx
import numpy as np

def load_graph(path):
    """
    Loads data from TSV file pointed to by path into a networkx graph
    """
    G = nx.Graph()
    with open(path, 'r') as inputFile:
        inputFile.readline()
        for line in inputFile:
            #if line.strip() == "#":
            #    break
            data = line.strip().split(" ")
            #print("new edge:", data[0], data[1], data[2])
            G.add_edge(data[0], data[1], confidence=data[2])
    return G



def load_disease_genes(path):
    """
    Loads disease genes from TSV file and returns a python list of all names
    """
    with open(path, 'r') as inputFile:
        return inputFile.read().splitlines()

def load_start_vector(path, PPI_Graph):
    diseaseGeneList = set(load_disease_genes(path))
    start_vector = np.zeros(PPI_Graph.number_of_nodes())
    numDiseaseGenes = len(diseaseGeneList)
    for i, node in enumerate(PPI_Graph.nodes()):
        if node in diseaseGeneList:
            start_vector[i] = 1/numDiseaseGenes
    return np.asarray(start_vector, order='F')
