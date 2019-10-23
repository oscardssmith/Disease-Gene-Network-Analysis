import networkx as nx
import matplotlib.pyplot as plt


def load_graph(path):
    G = nx.Graph()
    inputFile = open(path, 'r')
    for line in inputFile:
        if line.strip() == "#":
            break
        data = line.strip().split(" ")
        #print(data)
        for i in [0,1]:
            if data[i] not in nx.nodes(G):
                #print("new node:", data[i])
                G.add_node(data[i])
        #print("new edge:", data[0], data[1], data[2])
        G.add_edge(data[0], data[1], confidence=data[2])

    inputFile.close()

    return G








def main():
    PPI_Graph = load_graph("test-graph-data.tsv")

    nx.draw(PPI_Graph, node_color='r', edge_color='b')
    plt.show()


if __name__ == '__main__':
    main()
