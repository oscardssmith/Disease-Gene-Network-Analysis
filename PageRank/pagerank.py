import numpy as np
from scipy.spatial import distance

EPSILON = 0.001

# Given a np.array matrix, starting vector, prior bias vector, and back
# probability, calculate the rank of each node in the graph.
def pagerank(adjacency_matrix, starting_vector, prior_bias, beta):
    N = len(adjacency_matrix)
    sqrt_d_inverse = np.zeros((N, N))
    for i in range(N):
        sqrt_d_inverse[i][i] = 1/np.sqrt(sum(adjacency_matrix[i]))
    sqrt_d_inverse = np.array(sqrt_d_inverse)
    matrix = np.matmul(np.matmul(sqrt_d_inverse, adjacency_matrix), sqrt_d_inverse)

    d = float('inf')
    prev_vector = np.copy(starting_vector)
    iterations = 0
    while d > EPSILON:
        result = (1 - beta) * np.matmul(matrix, prev_vector)
        result = np.add(result, beta*prior_bias)
        d = distance.sqeuclidean(result, prev_vector)
        prev_vector = result
        iterations += 1
        print("iterations", iterations)
        print(result)
        print("difference", d)
    # zip prev_vector and gene names
    return prev_vector


def main():
    matrix = np.array([[0, 1, 1, 1], [1, 0, 1, 0], [1, 1, 0, 0], [1, 0, 0, 0]])
    starting_vector = np.array([1/4, 1/4, 1/4, 1/4])
    prior_bias = np.copy([0, 1/2, 0, 1/2])
    beta = 0.5
    result = pagerank(matrix, starting_vector, prior_bias, beta)
    print(result)

if __name__ == '__main__':
    main()
