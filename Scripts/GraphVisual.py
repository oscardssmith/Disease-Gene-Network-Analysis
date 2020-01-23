#import plotly.express as px
import matplotlib.pyplot as plt
import networkx as nx
from loader import load_graph

DATA_PATH = "../Data/9606.protein.links.v11.0.txt"

ppi_graph = load_graph(DATA_PATH)
result = []
hist = nx.degree_histogram(ppi_graph)
for i, count in enumerate(hist):
    if count > 0:
        result.append((i, count))
print(result)
plt.hist(hist, nbins=50)
plt.show()
