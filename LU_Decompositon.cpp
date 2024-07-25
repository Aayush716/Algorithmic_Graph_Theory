#include <iostream>
#include <vector>

using namespace std;

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
void printMatrix(const vector<vector<double>>& mat) {
    for (const auto& row : mat) {
        for (double val : row) {
            cout << val << "\t";
        }
        cout << endl;
    }
}

int main() {



    
    vector<vector<double>> A = {{2, -1, 0},
                                 {-1, 2, -1},
                                 {0, -1, 2}};
    
    vector<vector<double>> L, U;
    lu_decomposition(A, L, U);
    
    cout << "L:" << endl;
    printMatrix(L);
    cout << endl;
    
    cout << "U:" << endl;
    printMatrix(U);

    int mul1=1;
    int mul2=1;
    for(int i=0;i<L.size();i++){
        for(int j=0;j<L.size();j++){
            if(i==j){
                mul1 = mul1*L[i][j];
                mul2 = mul2*U[i][j];
            }
        }
    }
    cout<<"mul1 = "<<mul1<<" mul2 = "<<mul2;
    cout<<"matrix multiplicatin = "<<mul1*mul2;

    return 0;
}
