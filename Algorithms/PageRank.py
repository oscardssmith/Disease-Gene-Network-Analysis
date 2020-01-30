import sys
sys.path.insert(1, '../Scripts/')
#Insert relative paths for calls from run.py
sys.path.insert(1, 'Scripts/')
from CacheUtils import compute_if_not_cached
import GraphUtils
from CacheUtils import compute_if_not_cached
import os
from scipy.spatial import distance
import numpy as np
import math
from time import time
import networkx as nx
import loader
from loader import load_graph, load_start_vector

BETA = 0.3
DISEASE_GENE_FILE_PATH = "../Data/endometriosis-proteins.diseasegenes.tsv"
DATA_PATH = "../Data/9606.protein.links.v11.0.txt"
EPSILON = .000001 # 10^(-6)
PICKLE_PATH = "../Data/pickledmatrix"

def compute_matrix(graph):
    return np.asarray(GraphUtils.normalize_adjacency_matrix(
            nx.to_numpy_matrix(graph)))

# Given a np.array matrix, starting vector, prior bias vector, and back
# probability, calculate the rank of each node in the graph.
def rank_genes(graph, starting_vector, prior_bias, beta):
    #print("started ranking genes")
    start_time = time()

    # Load matrix from pickled object if exists to save time converting file.
    matrix = compute_if_not_cached(compute_matrix, graph)
    #print("time elapsed for normalizing the adjacency matrix: ", time() - start_time)

    d = float('inf')
    prev_vector = np.copy(starting_vector)
    iterations = 0
    while d > EPSILON:
        result = (1 - beta) * np.matmul(matrix, prev_vector)
        result = np.add(result, beta*prior_bias)
        d = distance.sqeuclidean(result, prev_vector)
        prev_vector = result
        iterations += 1
        #print("finished iterations", iterations)
    return prev_vector


def load_priors(priors_file, graph):
    prior_bias = np.zeros(graph.number_of_nodes())
    protein_list = []
    total = 0
    with open(priors_file, 'r') as inputFile:
        for line in inputFile:
            stripped = line.rstrip('\n')
            splitting = stripped.split('\t')
            if len(splitting) == 2:
                splitting[1] = int(splitting[1])
            protein_list.append(splitting)
    list_of_nodes = list(graph.nodes)
    for protein in protein_list:
        index = list_of_nodes.index(protein[0])
        if len(protein) == 1:
            prior_bias[index] = 1
            total += 1
        else:
            prior_bias[index] = protein[1]
            total += protein[1]
    prior_bias = (1/total) * prior_bias
    return prior_bias


def page_rank(graph, start_vector, prior_bias, beta=BETA):
    return GraphUtils.format_output(graph, rank_genes(graph, start_vector, prior_bias, beta))


def main():
    print(sys.argv)

    pathToPPINetworkFile = sys.argv[1]
    pathToDiseaseGeneFile = sys.argv[2]
    beta = float(sys.argv[3])

    print("loading data from files..")
    ppiGraph = compute_if_not_cached(load_graph, pathToPPINetworkFile, fileName="ppiGraph")
    diseaseGenes = load_start_vector(pathToDiseaseGeneFile, ppiGraph)

    return page_rank(ppiGraph, diseaseGenes, beta)

    '''print(time())
    # Read data from input file to networkx graph format.
    start_time = time()
    PPI_Graph = loader.load_graph(DATA_PATH)
    print("Graph loaded from file.\nTime elapsed:",
          time() - start_time, "seconds.")

    # Read data from disease gene file into list
    start_time = time()
    start_vector = loader.load_start_vector(DISEASE_GENE_FILE_PATH, PPI_Graph)
    print("Disease genes loaded from file.\nTime elapsed:",
          time() - start_time, "seconds.")

    start_time = time()
    probability_vector = PageRank(PPI_Graph, start_vector)
    print("PageRank matrix implementation finished running.\nTime elapsed:",
          time() - start_time, "seconds.")

    print(probability_vector)
    sorted_probabilities = sorted(probability_vector, key=lambda x: x[1])
    for i in range(50):
        print("gene name: ", sorted_probabilities[i][0],
              "probability:", sorted_probabilities[i][1])'''


if __name__ == '__main__':
    main()
