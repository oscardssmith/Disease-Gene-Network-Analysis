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
REAL_FILE = '../Data/9606.protein.links.v11.0.txt'

def diffusion_kernel(protein_graph, genes, B=1):
    print('in dk')
    L = nx.laplacian_matrix(protein_graph).todense()
    print('done densing')
    K = expm(-B * L)
    return np.matmul(K,  genes)

if __name__ == '__main__':
    t1 = time()
    #graph, label_to_idx = load_matrix(TEST_FILE, 0)
    graph = load_graph(REAL_FILE)
    print(time()-t1)
    #diffusion_kernel(graph, 1)
