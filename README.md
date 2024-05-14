# Algorithmic_Graph_Theory
This repo has solutions to standard graph problems

Generate a random graph on 11 vertices and find number of spanning trees in
it in a polynomial time

Pseudo code:

generateRandomGraph(n, p):
Initialize an empty adjacency matrix of size n x n
Seed the random number generator
Iterate over all pairs of nodes (i, j) with i < j:
Generate a random probability rand_prob between 0 and 1
If rand_prob is less than p:
Set adjacencyMatrix[i][j] = 1
Set adjacencyMatrix[j][i] = 1 (for undirected graph)
Return the adjacency matrix

laplacianMatrix(adjacencyMatrix):
Calculate the degree matrix:
Initialize a degree matrix of size n x n with zeros
Iterate over each node i:
Compute the degree of node i by summing its row in adjacencyMatrix
Set the diagonal element degreeMatrix[i][i] to the degree of node i

Compute the Laplacian matrix:
Initialize a Laplacian matrix of size n x n with zeros
Iterate over each node i and its adjacent node j:
Set Laplacian[i][j] = degreeMatrix[i][j] - adjacencyMatrix[i][j]
Return the Laplacian matrix
