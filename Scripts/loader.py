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



def load_test_graph():
    return nx.complete_graph(5)

def load_test_start_vector():
    return np.asarray([0,0,0,1,0])


def load_disease_genes(path):
    """
    Loads disease genes from TSV file and returns a python list of all names
    """
    with open(path, 'r') as inputFile:
        return inputFile.read().splitlines()

def load_start_vector(path, ppiGraph):
    diseaseGeneList = set(load_disease_genes(path))
    start_vector = np.zeros(ppiGraph.number_of_nodes())
    numDiseaseGenes = len(diseaseGeneList)
    for i, node in enumerate(ppiGraph.nodes()):
        if node in diseaseGeneList:
            start_vector[i] = 1/numDiseaseGenes
    return np.asarray(start_vector, order='F')
