import numpy as np
from scipy.linalg import expm

def diffusion_kernel(protein_adj_matrix, genes, B=1):
    A = protein_adj_matrix
    D = np.diag(A * np.ones(A.size()[1]))
    L = D-A
    K = expm(-B*L)
    return K*genes
    
def load_matrix(path)
    label_to_idx = {}
    with open(path) as edges:
        edges.
        
    with open(path) as edges:
    return graph, label_to_idx

if __name__ == '__main__':
    load_matrix('9606.protein.links.v11.0.tsv')
