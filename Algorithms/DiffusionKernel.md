# Diffusion Kernel

## Theory

The diffusion kernel method ranks genes by computing the matrix exponential of a laplacian matrix. For a PPI network, the laplacian matrix L is calculated as:

<img src="https://render.githubusercontent.com/render/math?math=L = D - A">

Where D is the diagonal matrix, and A is the adjacency matrix of the graph. We then obtain the diffusion kernel K by calculating the matrix exponential of L, multiplied by a magnitude of diffusion , a variable parameter:
<img src="https://render.githubusercontent.com/render/math?math=K = e^{L \times -\beta}">
To calculate the matrix exponential, we used eigenvector decomposition. Where:
L=PDP-1 
where P’s columns are the eigenvectors of L, and D is a diagonal matrix who’s entries are the eigenvalues of L. Since our graph is symmetric, the eigenvalues and eigenvectors are real, and can be efficiently calculated by BLAS. Furthermore, P is a normal matrix, so P-1=PT.Once we have done the decomposition, we can compute 
K = e(L  -)=Pe-βDPT
This is cheap because the exponential of a diagonal matrix is a diagonal matrix with it’s entries exponentiated. Further speedup of multiple runs can be obtained by storing and loading P and D to disk, which means that operations like finding optimal βcan be done very quickly.

With this method, we perform the same calculation, but with a much faster run time. To calculate the scores of disease genes, we multiply our starting gene vector, and the kernel, which results in a vector, where each element is a score that corresponds to a gene. This is done during the computation, rather than after, and thus the resulting vector is the desired list of scores.


![Random Walk with Restart Equation](../Images/RWR-equation.gif)


## Implementation

### symmetric_eigen_from_graph()

## diffusion_kernel_core()

### diffusion_kernel()

### main()
