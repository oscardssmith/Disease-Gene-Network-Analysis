# Diffusion Kernel

## Theory

The diffusion kernel method ranks genes by computing the matrix exponential of a laplacian matrix. For a PPI network, the laplacian matrix L is calculated as:

<img src="https://render.githubusercontent.com/render/math?math=L = D - A">

Where D is the diagonal matrix, and A is the adjacency matrix of the graph. We then obtain the diffusion kernel K by calculating the matrix exponential of L, multiplied by a magnitude of diffusion , a variable parameter:

<img src="https://render.githubusercontent.com/render/math?math=K = e^{(L \times -\beta)}">

To calculate the matrix exponential, we used eigenvector decomposition. Where:

<img src="https://render.githubusercontent.com/render/math?math=L = PDP^{-1}">

where P’s columns are the eigenvectors of L, and D is a diagonal matrix who’s entries are the eigenvalues of L. Since our graph is symmetric, the eigenvalues and eigenvectors are real, and can be efficiently calculated by BLAS. Furthermore, P is a normal matrix, so <img src="https://render.githubusercontent.com/render/math?math=P^{-1}=P^{T}">.Once we have done the decomposition, we can compute

<img src="https://render.githubusercontent.com/render/math?math=K = e^{(L \times -\beta)} = Pe^{\beta D}P^{T}">

This is cheap because the exponential of a diagonal matrix is a diagonal matrix with it’s entries exponentiated. Further speedup of multiple runs can be obtained by storing and loading P and D to disk, which means that operations like finding optimal βcan be done very quickly.

With this method, we perform the same calculation, but with a much faster run time. To calculate the scores of disease genes, we multiply our starting gene vector, and the kernel, which results in a vector, where each element is a score that corresponds to a gene. This is done during the computation, rather than after, and thus the resulting vector is the desired list of scores.

## Implementation

### symmetric_eigen_from_graph()
Computes the laplacian of the input PPI graph and returns its eigenvalues and eigenvectors.

### diffusion_kernel_core()
Performs the actual computation of the diffusion kernel and returns a sorted, formatted output vector containing each gene and its score.

### diffusion_kernel()
Takes three parameters:

- ppiGraph:
A networkx matrix representation of the PPI network to compute the diffusion kernel of.

- diseaseGenes:
A 1 x N vector (numpy array), where N is the number of genes in the PPI network. The value of the entry of each gene from the disease gene set is 1/number of genes in the disease gene set. Thus, each disease all share an equal probability.

- beta:
The magnitude of diffusion. Default value is 1.

### main()
The main method allows the wrapper method diffusion_kernel() to be run from the run.py in the command line. Command line arguments are parsed and passed to diffusion_kernel() and the output of diffusion_kernel() is written to a .csv file that is saved to the given file path.
