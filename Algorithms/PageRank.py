import sys
import os
sys.path.insert(1, '../Imports/')
# Insert relative paths for calls from run.py
sys.path.insert(1, 'Imports/')
import loader
import networkx as nx
from time import time
import math
import csv
import numpy as np
from scipy.spatial import distance
import GraphUtils
from CacheUtils import compute_if_not_cached

BETA = 0.3
DISEASE_GENE_FILE_PATH = "../Data/endometriosis-proteins.diseasegenes.tsv"
DATA_PATH = "../Data/9606.protein.links.v11.0.txt"
EPSILON = .000001  # 10^(-6)
PICKLE_PATH = "../Data/pickledmatrix"


def compute_matrix(graph):
    return np.asarray(GraphUtils.normalize_adjacency_matrix(
        nx.to_numpy_matrix(graph)))

# Given a np.array matrix, starting vector, prior bias vector, and back
# probability, calculate the rank of each node in the graph.


def rank_genes(graph, startingVector, priorBias, beta):
    print("Starting PageRank")

    # Load matrix from pickled object if exists to save time converting file.
    matrix = compute_if_not_cached(compute_matrix, graph)

    d = float('inf')
    prevVector = np.copy(startingVector)
    iterations = 0
    while d > EPSILON:
        result = (1 - beta) * np.matmul(matrix, prevVector)
        result = np.add(result, beta*priorBias)
        d = distance.sqeuclidean(result, prevVector)
        prevVector = result
        iterations += 1
        print("finished iteration:", iterations)
    return prevVector


def load_priors(priorsFile, graph):
    priorBias = np.zeros(graph.number_of_nodes())
    proteinList = []
    total = 0
    with open(priorsFile, 'r') as inputFile:
        for line in inputFile:
            stripped = line.rstrip('\n')
            splitting = stripped.split('\t')
            if len(splitting) == 2:
                splitting[1] = int(splitting[1])
            proteinList.append(splitting)
    listOfNodes = list(graph.nodes)
    for protein in proteinList:
        index = listOfNodes.index(protein[0])
        if len(protein) == 1:
            priorBias[index] = 1
            total += 1
        else:
            priorBias[index] = protein[1]
            total += protein[1]
    priorBias = (1/total) * priorBias
    return priorBias


def page_rank(graph, startVector, priorBias, beta=BETA):
    return GraphUtils.format_output(graph, rank_genes(graph, startVector, priorBias, beta))


def main():

    pathToPPINetworkFile = sys.argv[1]
    pathToDiseaseGeneFile = sys.argv[2]
    beta = float(sys.argv[3])
    outputFile = sys.argv[4]

    print("loading data from files..")
    ppiGraph = compute_if_not_cached(
        loader.load_graph, pathToPPINetworkFile, fileName="ppiGraph")
    diseaseGenes = loader.load_start_vector(pathToDiseaseGeneFile, ppiGraph)

    results = page_rank(ppiGraph, diseaseGenes, beta)

    print("Saving results to", outputFile)
    with open(outputFile, "w", newline='') as of:
        outputWriter = csv.writer(of, quoting=csv.QUOTE_ALL)
        outputWriter.writerow(["Gene", "Ranking"])
        for row in results:
            outputWriter.writerow(row)
    print("done.")

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
