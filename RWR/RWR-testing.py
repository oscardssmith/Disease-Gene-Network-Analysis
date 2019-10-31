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



if __name__ == '__main__':
    main()
