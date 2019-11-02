import numpy as np

# Calculates the euclidean distance between two vectors
def difference(v1, v2):
    # v1 list of nodes
    # v2 list of nodes
    magnitude = 0
    for i in range(len(v1)):
        magnitude += (v1[i] - v2[i])**2
    magnitude = np.sqrt(magnitude)
    return magnitude

# Given a np.array matrix, starting vector, prior bias vector, and back
# probability, calculate the rank of each node in the graph.
def pagerank(matrix, starting_vector, prior_bias, beta, epsilon):
    # return all rankings
    d = float('inf')
    prev_vector = np.copy(starting_vector)
    iterations = 0
    while d > epsilon:
        result = (1 - beta) * np.matmul(matrix, prev_vector)
        result = np.add(result, beta*prior_bias)
        d = difference(result, prev_vector)
        prev_vector = result
        iterations += 1
        print("iterations", iterations)
        print(result)
        print("difference", d)
    return prev_vector

def main():
    matrix = np.array([[0, 1/2, 1/2, 1], [1/3, 0, 1/2, 0], [1/3, 1/2, 0, 0], [1/3, 0, 0, 0]])
    starting_vector = np.array([1/4, 1/4, 1/4, 1/4])
    prior_bias = np.copy([0, 1/2, 0, 1/2])
    beta = 0.5
    epsilon = 0.001
    result = pagerank(matrix, starting_vector, prior_bias, beta, epsilon)
    print(result)

if __name__ == '__main__':
    main()
