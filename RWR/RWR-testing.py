import networkx as nx
import matplotlib.pyplot as plt
import time



def load_graph(path):
    """
    Loads data from TSV file pointed to by path into a networkx graph
    """
    G = nx.Graph()
    inputFile = open(path, 'r')
    for line in inputFile:
        #if line.strip() == "#":
        #    break
        data = line.strip().split(" ")
        #print("new edge:", data[0], data[1], data[2])
        G.add_edge(data[0], data[1], confidence=data[2])

    inputFile.close()

    return G










def main():

    #Read data from input file to networkx graph format.
    pathToData = "../test-graph-data.tsv"
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


if __name__ == '__main__':
    main()
