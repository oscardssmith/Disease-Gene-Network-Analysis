import sys
sys.path.insert(1, '../Scripts/')
from GraphUtils import normalize_adjacency_matrix
import loader
import networkx as nx
import matplotlib.pyplot as plt
import time
import math
import numpy as np
from scipy.spatial import distance


BETA = 0.3
EPSILON = 0.001

# Given a np.array matrix, starting vector, prior bias vector, and back
# probability, calculate the rank of each node in the graph.
def rank_genes(adjacency_matrix, starting_vector, prior_bias, beta):
    print("started ranking genes")
    startTime = time.time()
    matrix = normalize_adjacency_matrix(adjacency_matrix)
    endTime = time.time()
    print("time elapsed for normalizing the adjacency matrix: ",endTime - startTime )
    d = float('inf')
    prev_vector = np.copy(starting_vector)
    print("starting vector shape:", prev_vector.shape)
    iterations = 0
    while d > EPSILON:
        result = (1 - beta) * np.matmul(matrix, prev_vector)
        print("shape of result vector after matmul:", result.shape)
        result = np.add(result, beta*prior_bias)
        print("shape of result after np add:", result.shape)
        d = distance.sqeuclidean(result, prev_vector)
        prev_vector = result
        iterations += 1
        print("finished iterations", iterations)
        print("shape of the previous vector:", prev_vector.shape)
    return prev_vector

def PageRank(graph, start_vector):
    adjacency_matrix = nx.to_numpy_matrix(graph)
    prior_bias = np.copy(start_vector)
    print("adjacency matrix shape:", adjacency_matrix.shape)
    print("start vector shape:", start_vector.shape)
    return zip(graph.nodes(),rank_genes(adjacency_matrix, start_vector, prior_bias, BETA))

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
    for name, p in sortedProbabilities:
        print("gene name:", name, "probability", p)


if __name__ == '__main__':
    main()
