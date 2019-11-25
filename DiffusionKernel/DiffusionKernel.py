import sys
sys.path.insert(1, '../Scripts/')
from loader import load_graph, load_start_vector

import gzip
from time import time

import networkx as nx
import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

#import ..RWR.RandomWalk import load_graph

TEST_FILE = '../Data/fake.tsv'
REAL_FILE = '../Data/9606.protein.links.v11.0.txt'
pathToDiseaseGeneFile = "../Data/EndometriosisProteins.tsv"

def diffusion_kernel(PPI_Graph, genes, B=1):
    # np.asarray needed due to this bug https://github.com/scipy/scipy/issues/5546
    # Note that this bug is fixed, but not yet backported
    L = np.asarray(nx.laplacian_matrix(PPI_Graph, weight=-B).todense())
    print(L)
    return expm_multiply(L, genes)
    #K = expm(-B * L)
    #print(K.shape)
    #return np.matmul(K, genes)

if __name__ == '__main__':
    t1 = time()
    #PPI_Graph = load_graph(REAL_FILE)
    #nx.write_weighted_edgelist(PPI_Graph, 'graph.edgelist')
    PPI_Graph = nx.read_weighted_edgelist('graph.edgelist')
    start_vector = load_start_vector(pathToDiseaseGeneFile, PPI_Graph)
    t2 = time()
    print(t2-t1)
    x = diffusion_kernel(PPI_Graph, start_vector)
    print(x)
    print(time()-t2)
