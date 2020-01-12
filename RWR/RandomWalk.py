import sys
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

import pickle
import os






# Runs Random Walk with Restart using a matrix implementation
def randomWalkMatrix(matrix, start_vector, R, max_iterations, norm_threshold):

    print("STARTING RANDOM WALK")

    previous_vector = np.copy(start_vector)
    iterations = 0
    diff = float('inf')

    while diff > norm_threshold and iterations < max_iterations:
        print("iteration:", iterations)

        #Perform one step of the walk
        new_vector = (1 - R) * np.matmul(matrix, previous_vector)
        new_vector = np.add(new_vector, R * start_vector)


        diff = distance.sqeuclidean(new_vector, previous_vector)
        previous_vector = new_vector
        iterations += 1

    return new_vector








def RandomWalk(graph, start_vector):
    print("INITIALIZING RANDOM WALK")

    """
    This method can be called from anywhere (such as validation scripts) and does whatever it needs to do to produce a properly formatted output,
    using only the given parameters.

    @param graph: a networkx graph object containing the entire PPI network
    @param diseaseGeneList: a python list object where each item is a string containing the name of a known disease gene

    @returns: a nested list of tuples, in sorted order of probability, where each item contains the name of a gene, and its respective probability as determined by the algorithm
    """
    # Set algorithm constants
    R = 0.2
    max_iterations = 500
    norm_threshold = 10**(-6)
    print("creating matrix")

    # Load matrix from pickled object if exists to save time converting file.
    if os.path.isfile("pickledmatrix"):
        with open("pickledmatrix", 'r') as handle:
            matrix = pickle.load(handle)
    else:
        matrix = np.asarray(normalize_adjacency_matrix(nx.to_numpy_matrix(graph)))
        with open("pickledmatrix", 'w') as handle:
            pickle.dump(matrix, handle)


    probabilityVector = randomWalkMatrix(matrix, start_vector, R, max_iterations, norm_threshold)

    #format probabilityVector into usable output
    print("formatting output")
    return format_output(graph, probabilityVector)






def main():
    pathToData = "../Data/9606.protein.links.v11.0.txt"
    pathToDiseaseGeneFile = "../Data/EndometriosisProteins.tsv"


    print("Loading graph from file:", pathToData)

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
    probabilityVector = RandomWalk(PPI_Graph, start_vector)
    endTime = time.time()

    print("Random Walk matrix implementation finished running.\nTime elapsed:", endTime - startTime, "seconds.")
    print(probabilityVector, start_vector)


    #Visualize graph in matplotlib.
    #startTime = time.time()
    #nx.draw(PPI_Graph, node_color='r', edge_color='b')
    #endTime = time.time()

    #print("graph visualized.\nTime elapsed:", endTime - startTime, "seconds.")
    #plt.show()


    #Export to graphML file
    #startTime = time.time()
    #nx.write_graphml(PPI_Graph, "PPI_Network.graphml")
    #endTime = time.time()

    #print("graph exported.\nTime elapsed:", endTime - startTime, "seconds.")



if __name__ == '__main__':
    main()
