#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <algorithm>
using namespace std;

long long baseToDecimal(const string &value, int base) {
    long long result = 0, power = 1;
    for (int i = value.size() - 1; i >= 0; i--) {
        char c = value[i];
        int digit;
        if (isdigit(c)) digit = c - '0';
        else if (isalpha(c)) digit = (tolower(c) - 'a') + 10;
        else throw runtime_error("Invalid char in base value");
        if (digit >= base) throw runtime_error("Digit exceeds base");
        result += digit * power;
        power *= base;
    }
    return result;
}

// Gaussian elimination
vector<double> solveGaussian(vector<vector<double>> A, vector<double> b) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        int pivot = i;
        for (int j = i + 1; j < n; j++)
            if (fabs(A[j][i]) > fabs(A[pivot][i])) pivot = j;
        swap(A[i], A[pivot]);
        swap(b[i], b[pivot]);

        for (int j = i + 1; j < n; j++) {
            double factor = A[j][i] / A[i][i];
            for (int k = i; k < n; k++) A[j][k] -= factor * A[i][k];
            b[j] -= factor * b[i];
        }
    }

    vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        double sum = b[i];
        for (int j = i + 1; j < n; j++) sum -= A[i][j] * x[j];
        x[i] = sum / A[i][i];
    }
    return x;
}

// Solve polynomial and return constant term c
double solvePolynomial(const vector<pair<long long, long long>> &points, int k) {
    vector<vector<double>> A(k, vector<double>(k));
    vector<double> b(k);

    for (int i = 0; i < k; i++) {
        long long xi = points[i].first;
        long long yi = points[i].second;
        b[i] = yi;
        long long power = 1;
        for (int j = k - 1; j >= 0; j--) {
            A[i][j] = power;
            power *= xi;
        }
    }

    auto coeffs = solveGaussian(A, b);
    return coeffs.back(); // c (constant term)
}

int main() {
    cout << fixed << setprecision(2);

    // ---------------- Test Case 1 ----------------
    {
        int k = 3;
        vector<pair<long long,long long>> points;
        points.push_back({1, baseToDecimal("4", 10)});
        points.push_back({2, baseToDecimal("111", 2)});
        points.push_back({3, baseToDecimal("12", 10)});
        double c = solvePolynomial(points, k);
        cout << "Test Case 1 -> c = " << c << endl;
    }

    // ---------------- Test Case 2 ----------------
    {
        int k = 7;
        vector<pair<long long,long long>> points;
        points.push_back({1, baseToDecimal("13444211440455345511", 6)});
        points.push_back({2, baseToDecimal("aed7015a346d635", 15)});
        points.push_back({3, baseToDecimal("6aeeb69631c227c", 15)});
        points.push_back({4, baseToDecimal("e1b5e05623d881f", 16)});
        points.push_back({5, baseToDecimal("316034514573652620673", 8)});
        points.push_back({6, baseToDecimal("2122212201122002221120200210011020220200", 3)});
        points.push_back({7, baseToDecimal("20120221122211000100210021102001201112121", 3)});
        // we only need k points = 7, so stop here

        double c = solvePolynomial(points, k);
        cout << "Test Case 2 -> c = " << c << endl;
    }

    return 0;
}
