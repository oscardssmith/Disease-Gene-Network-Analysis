import sys
sys.path.insert(1, '../Imports/')
# Insert relative paths for calls from run.py
sys.path.insert(1, 'Imports/')
from CacheUtils import compute_if_not_cached
from GraphUtils import format_output
import networkx as nx
import numpy as np
import loader


def symmetric_eigen_from_graph(ppiGraph):
    L = nx.laplacian_matrix(ppiGraph).todense()
    return np.linalg.eigh(L)


def diffusion_kernel_core(ppiGraph, genes, beta):
    # Compute matrix exponential with eigen decomposition
    # Faster since it uses the fact that the matrix is real, symetric

    vals, vecs = compute_if_not_cached(symmetric_eigen_from_graph, ppiGraph)
    result = np.dot(np.dot(np.dot(genes, np.transpose(vecs)),
                           np.diag(np.exp(-beta*vals))), vecs)
    result = np.array(result).flatten()
    return format_output(ppiGraph, result)

# Standard Wrapper Function


def diffusion_kernel(ppiGraph, diseaseGenes, beta=1):
    print("running diffusion kernel..")
    return diffusion_kernel_core(ppiGraph, diseaseGenes, beta)


if __name__ == '__main__':

    pathToPPINetworkFile = sys.argv[1]
    pathToDiseaseGeneFile = sys.argv[2]
    beta = float(sys.argv[3])
    outputFile = sys.argv[4]

    print("loading data from files..")
    ppiGraph = compute_if_not_cached(
        loader.load_graph, pathToPPINetworkFile, fileName="ppiGraph")
    diseaseGenes = loader.load_start_vector(pathToDiseaseGeneFile, ppiGraph)

    results = diffusion_kernel(ppiGraph, diseaseGenes, beta)

    print("Saving results to", outputFile)
    with open(outputFile, "w", newline='') as of:
        outputWriter = csv.writer(of, quoting=csv.QUOTE_ALL)
        outputWriter.writerow(["Gene", "Ranking"])
        for row in results:
            outputWriter.writerow(row)
    print("done.")
