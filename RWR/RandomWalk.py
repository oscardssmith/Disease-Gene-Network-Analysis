
import sys
sys.path.insert(1, '../Scripts/')
import loader

import networkx as nx
import matplotlib.pyplot as plt
import time


def load_graph(path):
    """
    Loads data from TSV file pointed to by path into a networkx graph
    """
    G = nx.Graph()
    with open(path, 'r') as inputFile:
        inputFile.readline()
        for line in inputFile:
            #if line.strip() == "#":
            #    break
            data = line.strip().split(" ")
            #print("new edge:", data[0], data[1], data[2])
            G.add_edge(data[0], data[1], confidence=data[2])
    return G

def main():
    pathToData = "../test-graph-data.tsv"

    #Read data from input file to networkx graph format.
    startTime = time.time()
    PPI_Graph = load_graph(pathToData)
    endTime = time.time()

    print("graph loaded from file.\nTime elapsed:", endTime - startTime, "seconds.")


    #Visualize graph in matplotlib.
    startTime = time.time()
    nx.draw(PPI_Graph, node_color='r', edge_color='b')
    endTime = time.time()

    print("graph visualized.\nTime elapsed:", endTime - startTime, "seconds.")
    plt.show()


    #export to graphML file
    startTime = time.time()
    nx.write_graphml(PPI_Graph, "PPI_Network.graphml")
    endTime = time.time()

    print("graph exported.\nTime elapsed:", endTime - startTime, "seconds.")



def RandomWalk(graph, diseaseGeneList):
    """
    This method can be called from anywhere (such as validation scripts) and does whatever it needs to do to produce a properly formatted output,
    using only the given parameters.

    @param graph: a networkx graph object containing the entire PPI network
    @param diseaseGeneList: a python list object where each item is a string containing the name of a known disease gene

    @returns: a nested list of tuples, in sorted order of probability, where each item contains the name of a gene, and its respective probability as determined by the algorithm
    """

    return [['9606.ENSP00000351407', 0.45], ['genename4', 0.36], ['genename3', 0.301], ['genename2', 0.23]]



if __name__ == '__main__':
    main()
