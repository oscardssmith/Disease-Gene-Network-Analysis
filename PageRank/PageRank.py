import sys
sys.path.insert(1, '../Scripts/')
import GraphUtils
import os
from scipy.spatial import distance
import numpy as np
import math
import time
import matplotlib.pyplot as plt
import networkx as nx
import loader
try:
    import cPickle as pickle
except:
    import pickle


BETA = 0.3
DISEASE_GENE_FILE_PATH = "../Data/EndometriosisProteins.tsv"
DATA_PATH = "../Data/9606.protein.links.v11.0.txt"
EPSILON = 10**(-6)
PICKLE_PATH = "../Data/pickledmatrix"

# Given a np.array matrix, starting vector, prior bias vector, and back
# probability, calculate the rank of each node in the graph.


def rank_genes(adjacency_matrix, starting_vector, prior_bias, beta):
    print("started ranking genes")
    start_time = time.time()

    # Load matrix from pickled object if exists to save time converting file.
    if os.path.isfile(PICKLE_PATH):
        print("pickled matrix file exists, loading matrix from file")
        with open(PICKLE_PATH, 'rb') as handle:
            matrix = np.asarray(pickle.load(handle))
    else:
        matrix = np.asarray(GraphUtils.normalize_adjacency_matrix(
            nx.to_numpy_matrix(graph)))
        with open(PICKLE_PATH, 'wb') as handle:
            pickle.dump(matrix, handle)
    print("time elapsed for normalizing the adjacency matrix: ", time.time() - start_time)

    d = float('inf')
    prev_vector = np.copy(starting_vector)
    iterations = 0
    while d > EPSILON:
        result = (1 - beta) * np.matmul(matrix, prev_vector)
        result = np.add(result, beta*prior_bias)
        d = distance.sqeuclidean(result, prev_vector)
        prev_vector = result
        iterations += 1
        print("finished iterations", iterations)
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


def PageRank(graph, start_vector, prior_bias):
    adjacency_matrix = nx.to_numpy_matrix(graph)
    return GraphUtils.format_output(graph, rank_genes(adjacency_matrix, start_vector, prior_bias, BETA))


def main():

    # Read data from input file to networkx graph format.
    start_time = time.time()
    PPI_Graph = loader.load_graph(DATA_PATH)
    print("Graph loaded from file.\nTime elapsed:",
          time.time() - start_time, "seconds.")

    # Read data from disease gene file into list
    start_time = time.time()
    start_vector = loader.load_start_vector(DISEASE_GENE_FILE_PATH, PPI_Graph)
    print("Disease genes loaded from file.\nTime elapsed:",
          time.time() - start_time, "seconds.")

    start_time = time.time()
    probability_vector = PageRank(PPI_Graph, start_vector)
    print("PageRank matrix implementation finished running.\nTime elapsed:",
          time.time() - start_time, "seconds.")

    print(probability_vector)
    sorted_probabilities = sorted(probability_vector, key=lambda x: x[1])
    for i in range(50):
        print("gene name: ", sorted_probabilities[i][0],
              "probability:", sorted_probabilities[i][1])


if __name__ == '__main__':
    main()
