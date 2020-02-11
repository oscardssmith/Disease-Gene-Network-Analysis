import sys
import os
sys.path.insert(1, '../Imports/')
# Insert relative paths for calls from run.py
sys.path.insert(1, 'Imports/')
from scipy.spatial import distance
import numpy as np
import math
import time
import csv
import networkx as nx
import GraphUtils
from CacheUtils import compute_if_not_cached
import loader





def random_walk_matrix(matrix, startVector, R, maxIterations, normThreshold):
    """
    Runs Random Walk with Restart using a matrix implementation

    @param matrix: numpy array, normalized adjancency matrix of entire PPI network
    @param startVector: numpy array, contains weighted start probabilities
    @param R: float, probability of restart parameter
    @param maxInterations: integer, maximum number of iterations to run
    @param normThreshold: integer, threshold at which the algorithm stops running if the difference between two steps is less than it

    @returns numpy array, final vector containing ranked proteins
    """
    print("STARTING RANDOM WALK")

    previousVector = np.copy(startVector)
    iterations = 0
    diff = float('inf')

    while diff > normThreshold and iterations < maxIterations:
        print("iteration:", iterations)

        # Perform one step of the walk
        newVector = (1 - R) * np.matmul(matrix, previousVector)
        newVector = np.add(newVector, R * startVector)

        diff = distance.sqeuclidean(newVector, previousVector)
        previousVector = newVector
        iterations += 1

    return newVector


def create_normalized_matrix(ppiGraph):
    """
    Generates normalized adjacency matrix.

    @param ppiGraph: a networkx graph containing the entire PPI network
    @returns: a numpy array that contains the normalized adjacency matrix
    """

    return np.asarray(GraphUtils.normalize_adjacency_matrix(nx.to_numpy_matrix(ppiGraph)))


def random_walk(graph, startVector, r=0.4):

    """
    This method can be called from anywhere (such as validation scripts) and does whatever it needs to do to produce a properly formatted output,
    using only the given parameters.

    @param graph: a networkx graph object containing the entire PPI network
    @param startVector: a numpy array that contains the weighted start probabilities for each protein in the network

    @returns: a nested list of tuples, in sorted order of probability, where each item contains the name of a gene, and its respective probability as determined by the algorithm
    """

    print("INITIALIZING RANDOM WALK")

    maxIterations = 500
    normThreshold = 10**(-6)

    print("creating matrix")

    matrix = compute_if_not_cached(create_normalized_matrix, graph, fileName="rwr_normalized_matrix")

    probabilityVector = random_walk_matrix(matrix, startVector, r, maxIterations, normThreshold)

    # format probabilityVector into usable output
    print("formatting output")
    return GraphUtils.format_output(graph, probabilityVector)


def main():
    """
    Allows random_walk() to be run through run.py.
    Parses command line arguments and feeds them as parameters to random_walk().
    Outputs list of ranked proteins as a .csv file in specified file path. 
    """
    pathToPPINetworkFile = sys.argv[1]
    pathToDiseaseGeneFile = sys.argv[2]
    R = float(sys.argv[3])
    outputFile = sys.argv[4]

    print("loading data from files..")
    ppiGraph = compute_if_not_cached(loader.load_graph, pathToPPINetworkFile, fileName="ppiGraph")
    diseaseGenes = loader.load_start_vector(pathToDiseaseGeneFile, ppiGraph)

    results = random_walk(ppiGraph, diseaseGenes, R)

    print("Saving results to", outputFile)
    with open(outputFile, "w", newline='') as of:
        outputWriter = csv.writer(of, quoting=csv.QUOTE_ALL)
        outputWriter.writerow(["Gene", "Ranking"])
        for row in results:
            outputWriter.writerow(row)
    print("done.")


if __name__ == '__main__':
    main()
