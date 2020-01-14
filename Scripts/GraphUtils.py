import numpy as np

def normalize_adjacency_matrix(adjacency_matrix):
    #ravel there so that diag is nx1 instead of 1xn
    diag = np.ravel(1/np.sqrt(np.sum(adjacency_matrix, axis=1)))
    sqrt_d_inverse = np.diag(diag)
    return np.asarray(np.matmul(np.matmul(sqrt_d_inverse, adjacency_matrix), sqrt_d_inverse))



def format_output(graph, rawOutputVector):
    #format probabilityVector into usable output
    output = list(zip(graph.nodes(), rawOutputVector))
    output.sort(key=lambda tup: tup[1], reverse=True)

    return output
