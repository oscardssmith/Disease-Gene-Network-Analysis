import numpy as np

def normalize_adjacency_matrix(adjacency_matrix):
    #ravel there so that diag is nx1 instead of 1xn
    diag = np.ravel(1/np.sqrt(np.sum(adjacency_matrix, axis=1)))
    sqrt_d_inverse = np.diag(diag)
    return np.matmul(np.matmul(sqrt_d_inverse, adjacency_matrix), sqrt_d_inverse)
