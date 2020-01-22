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
    # Compute matrix exponential with eigen decomposition
    # Faster since it uses the fact that the matrix is real, symetric
    L = nx.laplacian_matrix(PPI_Graph).todense()
    genes = np.array(genes, order='C')
    vals, vecs = np.linalg.eigh(-B*L)
    K = np.dot(np.dot(vecs, np.diag(np.exp(vals))), np.transpose(vecs))
    return np.matmul(K, genes)

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
