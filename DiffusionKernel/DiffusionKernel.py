import sys
sys.path.insert(1, '../Scripts/')
#Insert relative paths for calls from run.py
sys.path.insert(1, 'Scripts/')
from CacheUtils import compute_if_not_cached
from loader import load_graph, load_start_vector
from GraphUtils import format_output

import gzip
import os
import pickle
from time import time

import networkx as nx
import numpy as np

# BETA = 1
# pathToDiseaseGeneFile = "../Data/endometriosis-proteins.diseasegenes.tsv"
# REAL_FILE = '../Data/9606.protein.links.v11.0.txt'
# TEST_FILE = '../Data/fake.tsv'

def symmetric_eigen_from_graph(ppiGraph):
    L = nx.laplacian_matrix(ppiGraph).todense()
    return np.linalg.eigh(L)

def diffusion_kernel_core(ppiGraph, genes, beta):
    # Compute matrix exponential with eigen decomposition
    # Faster since it uses the fact that the matrix is real, symetric
    
    vals, vecs = compute_if_not_cached(symmetric_eigen_from_graph, ppiGraph)
    result = np.dot(np.dot(np.dot(genes, np.transpose(vecs)), np.diag(np.exp(-beta*vals))), vecs)
    result = np.array(result).flatten()
    return format_output(ppiGraph, result)


# unnecessary??
# def dk_leaveOneOut(startVector):
#     try:
#         ppiGraph = nx.read_weighted_edgelist('graph.edgelist')
#     except:
#         ppiGraph = load_graph(REAL_FILE)
#         nx.write_weighted_edgelist(ppiGraph, 'graph.edgelist')
#     x = diffusion_kernel(ppiGraph, startVector)
#     return x



# Standard Wrapper Function
def diffusion_kernel(ppiGraph, diseaseGenes, beta = 1):
    print("running diffusion kernel..")
    return diffusion_kernel_core(ppiGraph, diseaseGenes, beta)
    # Format output if necessary for validation?





#unnecessary?
# def dk_test():
#     try:
#     	ppiGraph = nx.read_weighted_edgelist('graph.edgelist')
#     except:
#     	ppiGraph = load_graph(REAL_FILE)
#     	nx.write_weighted_edgelist(ppiGraph, 'graph.edgelist')
#     start_vector = load_start_vector(pathToDiseaseGeneFile, ppiGraph)
#     x = diffusion_kernel(ppiGraph, start_vector)
#     return x
    


if __name__ == '__main__':

    pathToPPINetworkFile = sys.argv[1]
    pathToDiseaseGeneFile = sys.argv[2]
    beta = float(sys.argv[3])

    print("loading data from files..")
    ppiGraph = compute_if_not_cached(load_graph, pathToPPINetworkFile, fileName="ppiGraph")
    diseaseGenes = load_start_vector(pathToDiseaseGeneFile, ppiGraph)

    diffusion_kernel(ppiGraph, diseaseGenes, beta)



    # tick = time()
    # test = dk_test()
    # time_elapsed = time() - tick
    # print('Finished. Time taken: ', time_elapsed)
    # for i in range(100):
    # 	print(test[i])
