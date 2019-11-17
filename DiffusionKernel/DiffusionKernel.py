import sys
sys.path.insert(1, '../Scripts/')
from loader import load_graph

import gzip
from time import time

import networkx as nx
import numpy as np
from scipy.linalg import expm

#import ..RWR.RandomWalk import load_graph

TEST_FILE = '../Data/fake.tsv'
REAL_FILE = '../9606.protein.links.v11.0.tsv'

def diffusion_kernel(protein_graph, genes, B=1):
    L = nx.laplacian_matrix(protein_graph)
    K = expm(-B * L)
    return K * genes

if __name__ == '__main__':
    t1 = time()
    #graph, label_to_idx = load_matrix(TEST_FILE, 0)
    graph = load_graph(TEST_FILE)
    print(time()-t1)
    #diffusion_kernel(graph, 1)
