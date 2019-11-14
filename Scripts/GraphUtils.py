import numpy as np

def normalize_adjacency_matrix(adjacency_matrix):
    sqrt_d_inverse = np.diag(1/np.sqrt(np.sum(adjacency_matrix, axis=1)))
    return np.matmul(np.matmul(sqrt_d_inverse, adjacency_matrix), sqrt_d_inverse)
