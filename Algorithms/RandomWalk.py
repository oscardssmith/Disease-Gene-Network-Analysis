import sys
import os
sys.path.insert(1, '../Scripts/')
# Insert relative paths for calls from run.py
sys.path.insert(1, 'Scripts/')
from scipy.spatial import distance
import numpy as np
import math
import time
import networkx as nx
import GraphUtils
from CacheUtils import compute_if_not_cached
import loader
try:
    import cPickle as pickle
except:
    import pickle

    


# Runs Random Walk with Restart using a matrix implementation
def random_walk_matrix(matrix, startVector, R, maxIterations, normThreshold):
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
    return np.asarray(GraphUtils.normalize_adjacency_matrix(nx.to_numpy_matrix(ppiGraph)))


def random_walk(graph, startVector, r=0.4):
    print("INITIALIZING RANDOM WALK")

    """
    This method can be called from anywhere (such as validation scripts) and does whatever it needs to do to produce a properly formatted output,
    using only the given parameters.

    @param graph: a networkx graph object containing the entire PPI network
    @param diseaseGeneList: a python list object where each item is a string containing the name of a known disease gene

    @returns: a nested list of tuples, in sorted order of probability, where each item contains the name of a gene, and its respective probability as determined by the algorithm
    """
    maxIterations = 500
    normThreshold = 10**(-6)
    print("creating matrix")

    matrix = compute_if_not_cached(create_normalized_matrix, graph, fileName="rwr_normalized_matrix")

    probabilityVector = random_walk_matrix(matrix, startVector, r, maxIterations, normThreshold)

    # format probabilityVector into usable output
    print("formatting output")
    return GraphUtils.format_output(graph, probabilityVector)


def main():

    pathToPPINetworkFile = sys.argv[1]
    pathToDiseaseGeneFile = sys.argv[2]
    R = float(sys.argv[3])

    print("loading data from files..")
    ppiGraph = compute_if_not_cached(loader.load_graph, pathToPPINetworkFile, fileName="ppiGraph")
    diseaseGenes = loader.load_start_vector(pathToDiseaseGeneFile, ppiGraph)

    random_walk(ppiGraph, diseaseGenes, R)


if __name__ == '__main__':
    main()
