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
    matrix = normalize_adjacency_matrix(adjacency_matrix)
    d = float('inf')
    prev_vector = np.copy(starting_vector)
    iterations = 0
    while d > EPSILON:
        result = (1 - beta) * np.matmul(matrix, prev_vector)
        result = np.add(result, beta*prior_bias)
        d = distance.sqeuclidean(result, prev_vector)
        prev_vector = result
        iterations += 1
        print("iterations", iterations)
        print(result)
        print("difference", d)
    # zip prev_vector and gene names
    return prev_vector

def PageRank(graph, disease_gene_list):
    adjacency_matrix = np.adjacency_matrix(graph)
    disease_gene_set = set(disease_gene_list)
    starting_vector = []
    num_disease_genes = len(disease_gene_list)
    for node in graph.nodes():
        if node in disease_gene_set:
            starting_vector.append(1/num_disease_genes)
        else:
            starting_vector.append(0)
    starting_vector = np.array(starting_vector, order='F')
    prior_bias = np.copy(starting_vector)
    return zip(graph.nodes(),rank_genes(adjacency_matrix, starting_vector, prior_bias, BETA))

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
    diseaseGeneList = loader.load_disease_genes(pathToDiseaseGeneFile)
    endTime = time.time()

    print("Disease genes loaded from file.\nTime elapsed:", endTime - startTime, "seconds.")


    startTime = time.time()
    probabilityVector = PageRank(PPI_Graph, diseaseGeneList)
    endTime = time.time()

    print("PageRank matrix implementation finished running.\nTime elapsed:", endTime - startTime, "seconds.")
    print(probabilityVector)


if __name__ == '__main__':
    main()
