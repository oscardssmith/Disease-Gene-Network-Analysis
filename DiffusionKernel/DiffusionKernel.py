import sys
sys.path.insert(1, '../Scripts/')
from loader import load_graph, load_start_vector
from GraphUtils import format_output
import gzip
import os
import pickle
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
    
    if os.path.isfile("../Data/pickled_eigen_decomp"):
        print("pickled matrix file exists, loading matrix from file")
        with open("../Data/pickled_eigen_decomp", 'rb') as handle:
            vals, vecs = pickle.load(handle)
    else:
        L = nx.laplacian_matrix(PPI_Graph).todense()
        vals, vecs = np.linalg.eigh(L)
        with open("../Data/pickled_eigen_decomp", 'wb') as handle:
            pickle.dump((vals, vecs), handle)
    result = np.dot(np.dot(np.dot(genes, np.transpose(vecs)), np.diag(np.exp(-B*vals))), vecs)
    return np.array(result).flatten()

def dk_test():
    try:
    	PPI_Graph = nx.read_weighted_edgelist('graph.edgelist')
    except:
    	PPI_Graph = load_graph(REAL_FILE)
    	nx.write_weighted_edgelist(PPI_Graph, 'graph.edgelist')
    start_vector = load_start_vector(pathToDiseaseGeneFile, PPI_Graph)
    x = diffusion_kernel(PPI_Graph, start_vector)
    formatted_vector = format_output(PPI_Graph, x)
    return formatted_vector

def dk_leaveOneOut(startVector):
    try:
        PPI_Graph = nx.read_weighted_edgelist('graph.edgelist')
    except:
        PPI_Graph = load_graph(REAL_FILE)
        nx.write_weighted_edgelist(PPI_Graph, 'graph.edgelist')
    x = diffusion_kernel(PPI_Graph, start_vector)
    formatted_vector = format_output(PPI_Graph, x)
    return formatted_vector


if __name__ == '__main__':
    tick = time()
    test = dk_test()
    time_elapsed = time() - tick
    print('Finished. Time taken: ', time_elapsed)
    for i in range(100):
    	print(test[i])
