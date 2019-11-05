
import sys
sys.path.insert(1, '../Scripts/')
import loader

import networkx as nx
import matplotlib.pyplot as plt
import time
import numpy as np




# Calculates the euclidean distance between two vectors
# Same function as in PageRank algorithm
def difference(v1, v2):
    # v1 list of nodes
    # v2 list of nodes
    magnitude = 0
    for i in range(len(v1)):
        magnitude += (v1[i] - v2[i])**2
    magnitude = np.sqrt(magnitude)
    return magnitude


# Runs Random Walk with Restart using a matrix implementation
def randomWalkMatrix(matrix, start_vector, R, max_iterations, norm_threshold):

    previous_vector = np.copy(start_vector)
    iterations = 0
    diff = float('inf')

    while diff > norm_threshold and iterations < max_iterations:
        #Perform one step of the walk
        new_vector = (1 - R) * np.matmul(matrix, previous_vector)
        new_vector = np.add(new_vector, R * start_vector)

        print(new_vector)
        print(previous_vector)


        diff = difference(new_vector, previous_vector)
        previous_vector = new_vector
        iterations += 1

    return new_vector








def RandomWalk(graph, diseaseGeneList):
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
    matrix = np.array(nx.to_numpy_matrix(graph))

    #compute start vector from disease gene list
    start_vector = []
    numDiseaseGenes = len(diseaseGeneList)
    for node in graph.nodes():
        if node in diseaseGeneList:
            start_vector.append(1/numDiseaseGenes)
        else:
            start_vector.append(0)
    start_vector = np.array(start_vector)


    probabilityVector = randomWalkMatrix(matrix, start_vector, R, max_iterations, norm_threshold)
    print(probabilityVector)

    return probabilityVector





def main():
    pathToData = "../Data/test-graph-data.tsv"
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
    probabilityVector = RandomWalk(PPI_Graph, diseaseGeneList)
    endTime = time.time()

    print("Random Walk matrix implementation finished running.\nTime elapsed:", endTime - startTime, "seconds.")
    print(probabilityVector)


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
