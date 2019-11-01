import networkx as nx




def load_graph(path):
    """
    Loads data from TSV file pointed to by path into a networkx graph

    @param path: the file location of the TSV file
    @returns: a networkx graph built from the contents of the TSV file
    """
    G = nx.Graph()
    inputFile = open(path, 'r')
    for line in inputFile:
        data = line.strip().split(" ")
        G.add_edge(data[0], data[1], confidence=data[2])

    inputFile.close()

    return G
