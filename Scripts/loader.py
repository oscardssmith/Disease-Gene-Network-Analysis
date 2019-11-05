import networkx as nx

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
