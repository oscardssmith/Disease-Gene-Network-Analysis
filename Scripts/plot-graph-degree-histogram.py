#import plotly.express as px
import matplotlib.pyplot as plt
import networkx as nx
from loader import load_graph
import sys
import os

if len(sys.argv) != 2:
    print("Usage: python3 plot-graph-degree-histogram.py path-to-ppi-network")
    sys.exit()

DATA_PATH = sys.argv[1]

ppi_graph = load_graph(DATA_PATH)
result = []
hist = nx.degree_histogram(ppi_graph)
for i, count in enumerate(hist):
    if count > 0:
        result.append((i, count))
print(result)
plt.hist(hist, nbins=50)
plt.show()
