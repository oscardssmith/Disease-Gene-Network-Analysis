import sys
sys.path.insert(1, '../Scripts/')
from GraphUtils import normalize_adjacency_matrix, format_output
import loader
import networkx as nx
import matplotlib.pyplot as plt
import time
import math
import numpy as np
from scipy.spatial import distance
try:
   import cPickle as pickle
except:
   import pickle

import os


BETA = 0.3
EPSILON = 10**(-6)

# Given a np.array matrix, starting vector, prior bias vector, and back
# probability, calculate the rank of each node in the graph.
def rank_genes(adjacency_matrix, starting_vector, prior_bias, beta):
    print("started ranking genes")
    startTime = time.time()

    # Load matrix from pickled object if exists to save time converting file.
    if os.path.isfile("../Data/pickledmatrix"):
        print("pickled matrix file exists, loading matrix from file")
        with open("../Data/pickledmatrix", 'rb') as handle:
            matrix = np.asarray(pickle.load(handle))
    else:
        matrix = np.asarray(normalize_adjacency_matrix(nx.to_numpy_matrix(graph)))
        with open("../Data/pickledmatrix", 'wb') as handle:
            pickle.dump(matrix, handle)
    endTime = time.time()
    print("time elapsed for normalizing the adjacency matrix: ",endTime - startTime )


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
        protein_file = inputFile.read()
        for line in protein_file:
            protein_list.append(line.split('\t'))
    list_of_nodes = list(graph.nodes)
    for protein in protein_list:
        index = list_of_nodes.index(protein[0])
        if protein.length == 1:
            prior_bias[index] = 1
            total += 1
        else:
            prior_bias[index] = protein[1]
            total += protein[1]
    prior_bias = (1/total) * prior_bias
    return prior_bias

def PageRank(graph, start_vector, prior_bias):
    adjacency_matrix = nx.to_numpy_matrix(graph)
    return format_output(graph,rank_genes(adjacency_matrix, start_vector, prior_bias, BETA))

def main():
    pathToData = "../Data/9606.protein.links.v11.0.txt"
    pathToDiseaseGeneFile = "../Data/EndometriosisProteins.tsv"


    #Read data from input file to networkx graph format.
    startTime = time.time()
    PPI_Graph = loader.load_graph(pathToData)
    endTime = time.time()

    print("Graph loaded from file.\nTime elapsed:", endTime - startTime, "seconds.")


    #Read data from disease gene file into list
    startTime = time.time()
    start_vector = loader.load_start_vector(pathToDiseaseGeneFile, PPI_Graph)
    endTime = time.time()

    print("Disease genes loaded from file.\nTime elapsed:", endTime - startTime, "seconds.")


    startTime = time.time()
    probabilityVector = PageRank(PPI_Graph, start_vector)
    endTime = time.time()

    print("PageRank matrix implementation finished running.\nTime elapsed:", endTime - startTime, "seconds.")
    print(probabilityVector)
    sortedProbabilities = sorted(probabilityVector, key=lambda x: x[1])
    for i in range(50):
        print("gene name: ", sortedProbabilities[i][0], "probability:", sortedProbabilities[i][1])


if __name__ == '__main__':
    main()
