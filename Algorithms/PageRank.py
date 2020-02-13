import sys
import os
sys.path.insert(1, '../Imports/')
# Insert relative paths for calls from run.py
sys.path.insert(1, 'Imports/')
import loader
import networkx as nx
import csv
import numpy as np
from scipy.spatial import distance
import GraphUtils
from CacheUtils import compute_if_not_cached

BETA = 0.4
EPSILON = .000001  # 10^(-6)

def compute_matrix(graph):
    return np.asarray(GraphUtils.normalize_adjacency_matrix(
        nx.to_numpy_matrix(graph)))

# Given a np.array matrix, starting vector, prior bias vector, and back
# probability, calculate the rank of each node in the graph.
def rank_genes(graph, startingVector, priorBias, beta):
    print("Starting PageRank")

    # Load matrix from pickled object if exists to save time converting file.
    matrix = compute_if_not_cached(compute_matrix, graph, fileName=graph.name)

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
    print("Finished PageRank")
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


if __name__ == '__main__':
    main()
