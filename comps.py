import gzip
from time import time

import numpy as np
from scipy.linalg import expm

TEST_FILE = 'fake.tsv'

def diffusion_kernel(protein_adj_matrix, genes, B=1):
    A = protein_adj_matrix
    D = np.diag(A * np.ones(A.shape[1]))
    L = D-A
    K = expm(-B * L)
    return K * genes

def load_matrix(path, min_weight):
    label_to_idx = {}
    max_ind = 0
    links = []
    with open(path, 'rt') as edges:
        edges.readline()
        for line in edges:
            protein1, protein2, weight = line.split()
            for protein in (protein1, protein2):
                if not protein in label_to_idx:
                    label_to_idx[protein] = max_ind
                    max_ind += 1
            if int(weight) >= min_weight:
                links.append((label_to_idx[protein1], label_to_idx[protein2]))

    graph = np.zeros((max_ind, max_ind))
    for x, y in links:
        graph[x, y] = 1
        graph[y, x] = 1
    return graph, label_to_idx

if __name__ == '__main__':
    t1 = time()
    graph, label_to_idx = load_matrix(TEST_FILE, 0)
    print(graph)
    print(label_to_idx)
    print(graph.size)
    print(time()-t1)
