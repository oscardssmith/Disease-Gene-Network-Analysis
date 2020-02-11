import numpy as np
import StringNameConverter as snc


def normalize_adjacency_matrix(adjacency_matrix):
    # ravel there so that diag is nx1 instead of 1xn
    diag = np.ravel(1/np.sqrt(np.sum(adjacency_matrix, axis=1)))
    sqrt_d_inverse = np.diag(diag)
    return np.asarray(np.matmul(np.matmul(sqrt_d_inverse, adjacency_matrix), sqrt_d_inverse))


def format_output(graph, raw_output_vector):
    # format probabilityVector into usable output
    output = list(zip(graph.nodes(), raw_output_vector))
    output.sort(key=lambda tup: tup[1], reverse=True)
    print(output[0])
    table = snc.load_lookup_table()
    for tup in output:
        tup[0] = snc.string_to_name(table, tup[0])
    print(tup[0])
    print(tup[1])
    return output
