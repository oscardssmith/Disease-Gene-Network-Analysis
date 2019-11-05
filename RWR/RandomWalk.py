
import sys
sys.path.insert(1, '../Scripts/')
import loader

import networkx as nx
import matplotlib.pyplot as plt
import time
import numpy as np





def randomWalkMatrix(matrix, start_vector, R, max_iterations, norm_threshold):

    previous_vector = np.copy(start_vector)
    iterations = 0
    difference = float('inf')

    while difference > norm_threshold and iterations < max_iterations:
        #Perform one step of the walk
        new_vector = (1 - R) * np.matmul(matrix, previous_vector)
        new_vector = np.add(new_vector, R * start_vector)

        difference = np.difference(new_vector, previous_vector)
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


    return [['9606.ENSP00000351407', 0.45], ['genename4', 0.36], ['genename3', 0.301], ['genename2', 0.23]]





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
    diseaseGeneList = load_disease_genes(pathToDiseaseGeneFile)
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


    print("graph exported.\nTime elapsed:", endTime - startTime, "seconds.")



if __name__ == '__main__':
    main()
