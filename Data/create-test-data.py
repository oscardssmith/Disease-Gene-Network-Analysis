"""
This script was used to create the test data file by taking the STRING data for a specific organism and eliminating a few of the proteins to cut down the size of the network.
"""


import networkx as nx
import matplotlib.pyplot as plt


def load_graph(path):
    G = nx.Graph()
    inputFile = open(path, 'r')
    outputFile = open("test-graph-data.tsv", 'w')
    for line in inputFile:
        if line.strip() == "#":
            break
        data = line.strip().split(" ")
        if not (data[0] == "568816.Acin_0001" or data[1] == "568816.Acin_0001" or
        data[0] == "568816.Acin_0002" or data[1] == "568816.Acin_0002" or
        data[0] == "568816.Acin_0003" or data[1] == "568816.Acin_0003" or
        data[0] == "568816.Acin_0004" or data[1] == "568816.Acin_0004" or
        data[0] == "568816.Acin_0005" or data[1] == "568816.Acin_0005" or
        data[0] == "568816.Acin_0006" or data[1] == "568816.Acin_0006" or
        data[0] == "568816.Acin_0007" or data[1] == "568816.Acin_0007" or
        data[0] == "568816.Acin_0008" or data[1] == "568816.Acin_0008" or
        data[0] == "568816.Acin_0009" or data[1] == "568816.Acin_0009" or
        data[0] == "568816.Acin_0010" or data[1] == "568816.Acin_0010"):
            outputFile.write(line)

        #print(data)
        for i in [0,1]:
            if data[i] not in nx.nodes(G):
                #print("new node:", data[i])
                G.add_node(data[i])
        #print("new edge:", data[0], data[1], data[2])
        G.add_edge(data[0], data[1], confidence=data[2])

    inputFile.close()
    outputFile.close()

    return G




def main():
    PPI_Graph = load_graph("test-graph-data-old.tsv")

    nx.draw(PPI_Graph, node_color='r', edge_color='b')
    plt.show()


if __name__ == '__main__':
    main()
