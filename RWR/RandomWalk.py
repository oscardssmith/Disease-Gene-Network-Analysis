import sys
import os
sys.path.insert(1, '../Scripts/')
import loader
from GraphUtils import normalize_adjacency_matrix
from GraphUtils import format_output
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








# Runs Random Walk with Restart using a matrix implementation
def random_walk_matrix(matrix, startVector, R, maxIterations, normThreshold):

    #print("STARTING RANDOM WALK")

    previousVector = np.copy(startVector)
    iterations = 0
    diff = float('inf')

    while diff > normThreshold and iterations < maxIterations:
        #print("iteration:", iterations)

        #Perform one step of the walk
        newVector = (1 - R) * np.matmul(matrix, previousVector)
        newVector = np.add(newVector, R * startVector)


        diff = distance.sqeuclidean(newVector, previousVector)
        previousVector = newVector
        iterations += 1

    return newVector








def random_walk(graph, startVector, r=0.3):
    print("INITIALIZING RANDOM WALK")

    """
    This method can be called from anywhere (such as validation scripts) and does whatever it needs to do to produce a properly formatted output,
    using only the given parameters.

    @param graph: a networkx graph object containing the entire PPI network
    @param diseaseGeneList: a python list object where each item is a string containing the name of a known disease gene

    @returns: a nested list of tuples, in sorted order of probability, where each item contains the name of a gene, and its respective probability as determined by the algorithm
    """
    # Set algorithm constants-- R is set to be same as Beta in PageRank for comparison
    #R = 0.3
    maxIterations = 500
    normThreshold = 10**(-6)
    #print("creating matrix")

    # Load matrix from pickled object if exists to save time converting file.
    if os.path.isfile("../Data/pickledmatrix"):
        print("pickled matrix file exists, loading matrix from file")
        with open("../Data/pickledmatrix", 'rb') as handle:
            matrix = np.asarray(pickle.load(handle))
    else:
        matrix = np.asarray(normalize_adjacency_matrix(nx.to_numpy_matrix(graph)))
        with open("../Data/pickledmatrix", 'wb') as handle:
            pickle.dump(matrix, handle)


    probabilityVector = random_walk_matrix(matrix, startVector, r, maxIterations, normThreshold)

    #format probabilityVector into usable output
    #print("formatting output")
    return format_output(graph, probabilityVector)






def main():
    pathToData = "../Data/9606.protein.links.v11.0.txt"
    pathToDiseaseGeneFile = "../Data/EndometriosisProteins.tsv"


    print("Loading graph from file:", pathToData)

    #Read data from input file to networkx graph format.
    startTime = time.time()
    ppiGraph = loader.load_graph(pathToData)
    endTime = time.time()

    print("Graph loaded from file.\nTime elapsed:", endTime - startTime, "seconds.")


    #Read data from disease gene file into list
    startTime = time.time()
    startVector = loader.load_start_vector(pathToDiseaseGeneFile, ppiGraph)
    endTime = time.time()

    print("Disease genes loaded from file.\nTime elapsed:", endTime - startTime, "seconds.")


    startTime = time.time()
    probabilityVector = random_walk(ppiGraph, startVector)
    endTime = time.time()

    print("Random Walk matrix implementation finished running.\nTime elapsed:", endTime - startTime, "seconds.")
    print(probabilityVector, startVector)


    #Visualize graph in matplotlib.
    #startTime = time.time()
    #nx.draw(ppiGraph, node_color='r', edge_color='b')
    #endTime = time.time()

    #print("graph visualized.\nTime elapsed:", endTime - startTime, "seconds.")
    #plt.show()


    #Export to graphML file
    #startTime = time.time()
    #nx.write_graphml(ppiGraph, "PPI_Network.graphml")
    #endTime = time.time()

    #print("graph exported.\nTime elapsed:", endTime - startTime, "seconds.")



if __name__ == '__main__':
    main()
