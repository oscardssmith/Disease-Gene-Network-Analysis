import gzip
from time import time

import networkx as nx
import numpy as np
from scipy.linalg import expm

#import ..RWR.RandomWalk import load_graph

TEST_FILE = '../Data/fake.tsv'

def diffusion_kernel(protein_graph, genes, B=1):
    L = nx.laplacian_matrix(protein_graph)
    K = expm(-B * L)
    return K * genes


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

if __name__ == '__main__':
    t1 = time()
    #graph, label_to_idx = load_matrix(TEST_FILE, 0)
    graph = load_graph(TEST_FILE)
    diffusion_kernel(graph, 1)
    print(time()-t1)
