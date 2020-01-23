import sys
sys.path.insert(1, '../Scripts/')
from scipy.sparse.linalg import expm_multiply
from scipy.linalg import expm
import numpy as np
import networkx as nx
from time import time
import gzip
from loader import load_graph, load_start_vector

BETA = 1
DISEASE_GENE_FILE_PATH = "../Data/EndometriosisProteins.tsv"
REAL_FILE = '../Data/9606.protein.links.v11.0.txt'
TEST_FILE = '../Data/fake.tsv'

def diffusion_kernel(PPI_Graph, genes, beta=BETA):
    # Compute matrix exponential with eigen decomposition
    # Faster since it uses the fact that the matrix is real, symetric
    L = nx.laplacian_matrix(PPI_Graph).todense()
    genes = np.array(genes, order='C')
    vals, vecs = np.linalg.eigh(-beta*L)
    K = np.dot(np.dot(vecs, np.diag(np.exp(vals))), np.transpose(vecs))
    return np.matmul(K, genes)

if __name__ == '__main__':
    start_time = time()
    #PPI_Graph = load_graph(REAL_FILE)
    #nx.write_weighted_edgelist(PPI_Graph, 'graph.edgelist')
    PPI_Graph = nx.read_weighted_edgelist('graph.edgelist')
    start_vector = load_start_vector(DISEASE_GENE_FILE_PATH, PPI_Graph)
    t1 = time()
    print(t2-start_time)
    result = diffusion_kernel(PPI_Graph, start_vector)
    print(result)
    print(time()-t1)
