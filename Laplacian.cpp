#include <iostream>
#include <vector>
#include <bits/stdc++.h>
//#include "Random_graph.h"

using namespace std;


vector<vector<int>> generateRandomGraph(int n, double p) {
    vector<vector<int>> adjacencyMatrix(n, vector<int>(n, 0));

    // Seed the random number generator
    srand(time(nullptr));

    // Generate edges randomly based on probability p
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double rand_prob = static_cast<double>(rand()) / RAND_MAX;
            if (rand_prob < p) {
                adjacencyMatrix[i][j] = 1;
                adjacencyMatrix[j][i] = 1; // Undirected graph, so the adjacency matrix is symmetric
            }
        }
    }

    return adjacencyMatrix;
}

// Function to print the adjacency matrix
void printAdjacencyMatrix(const vector<vector<int>>& adjacencyMatrix) {
    for (const auto& row : adjacencyMatrix) {
        for (int val : row) {
            cout << val << " ";
        }
        cout << endl;
    }
}


// Function to compute the Laplacian matrix
vector<vector<double>> laplacianMatrix(const vector<vector<int>>& adjacencyMatrix) {
    int n = adjacencyMatrix.size();

    // Initialize degree matrix with zeroes
    vector<vector<int>> degreeMatrix(n, vector<int>(n, 0));

    // Compute degree matrix
    for (int i = 0; i < n; ++i) {
        int degree = 0;
        for (int j = 0; j < n; ++j) {
            degree += adjacencyMatrix[i][j];
        }
        degreeMatrix[i][i] = degree;
    }

    // Compute Laplacian matrix: L = D - A
    vector<vector<double>> laplacian(n, vector<double>(n, 0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            laplacian[i][j] = degreeMatrix[i][j] - adjacencyMatrix[i][j];
        }
    }

    return laplacian;
}


void lu_decomposition(vector<vector<double>>& matrix, vector<vector<double>>& L, vector<vector<double>>& U) {
    int n = matrix.size();
    L.resize(n, vector<double>(n, 0));
    U.resize(n, vector<double>(n, 0));

    for (int i = 0; i < n; ++i) {
        L[i][i] = 1;

        for (int j = 0; j <= i; ++j) {
            double sum_upper = 0;
            for (int k = 0; k < j; ++k)
                sum_upper += L[i][k] * U[k][j];
            U[i][j] = matrix[i][j] - sum_upper;
        }

        for (int j = i; j < n; ++j) {
            double sum_lower = 0;
            for (int k = 0; k < i; ++k)
                sum_lower += L[j][k] * U[k][i];
            L[j][i] = (matrix[j][i] - sum_lower) / U[i][i];
        }
    }
}

// Utility function to print matrix
void printMatrix(const vector<vector<double>>& matrix) {
    for (const auto& row : matrix) {
        for (int val : row) {
            cout << val << " ";
        }
        cout << endl;
    }
}

double determinant(const vector<vector<double>>& matrix) {
    int n = matrix.size();

    // Base case: if matrix is 1x1, return its only element
    if (n == 1) {
        return matrix[0][0];
    }

    double det = 0;

    // Iterate over the first row to compute the cofactor expansion
    for (int j = 0; j < n; ++j) {
        // Compute the cofactor for the current element
        vector<vector<double>> submatrix(n - 1, vector<double>(n - 1, 0));
        for (int row = 1; row < n; ++row) {
            int colIndex = 0;
            for (int col = 0; col < n; ++col) {
                if (col != j) {
                    submatrix[row - 1][colIndex++] = matrix[row][col];
                }
            }
        }

        // Compute the sign multiplier
        double sign = (j % 2 == 0) ? 1 : -1;

        // Recursive call to compute the determinant of the submatrix
        det += sign * matrix[0][j] * determinant(submatrix);
    }

    return det;
}




int main() {
    // Example adjacency matrix
    // vector<vector<int>> adjacencyMatrix = {
    //     {0, 1, 0, 0},
    //     {1, 0, 1, 1},
    //     {0, 1, 0, 1},
    //     {0, 1, 1, 0}
    // };

    // Compute Laplacian matrix
    cout<<"enter number of nodes"<<endl;
    int numNodes;
    cin>>numNodes;

    // Probability of having an edge between any pair of nodes
    double probability = 0.7;
    vector<vector<int>> adjacencyMatrix = generateRandomGraph(numNodes, probability);
    cout << "Adjacency matrix of the random graph:" << endl;
    printAdjacencyMatrix(adjacencyMatrix);
    vector<vector<double>> laplacian = laplacianMatrix(adjacencyMatrix);

    // Print Laplacian matrix
    cout << "Laplacian matrix:" << endl;
    printMatrix(laplacian);

    vector<vector<double>> cofactor_laplacian;
    for(int i=0;i<laplacian.size()-1;i++){
        vector<double> temp;
        for(int j=0;j<laplacian[i].size()-1;j++){
            temp.push_back(laplacian[i][j]);
        }
        cofactor_laplacian.push_back(temp);
    }
    cout<<endl<<endl;
    
    cout <<"Cofactor matrix of the Laplacian matrix:" << endl;
    for(int i=0;i<laplacian.size()-1;i++){
        for(int j=0;j<laplacian[i].size()-1;j++){
            cout<<cofactor_laplacian[i][j]<<"  ";
        }
        cout<<endl;
    }


    double det = determinant(cofactor_laplacian);

    cout << "Determinant: " << det << endl;

    return 0;
}
