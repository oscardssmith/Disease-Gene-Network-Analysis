#import plotly.express as px
import matplotlib.pyplot as plt
import networkx as nx
from loader import load_graph

pathToData = "../Data/9606.protein.links.v11.0.txt"
PPI_Graph = load_graph(pathToData)

result = []
hist = nx.degree_histogram(PPI_Graph)
for i, count in enumerate(hist):
    if count>0:
        result.append((i,count))
print(result)
plt.hist(hist, nbins=50)
plt.show()
