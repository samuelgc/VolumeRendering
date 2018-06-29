#include "Matrix.h"


Matrix::Matrix() {
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            values[i][j] = 0;
        }
    }
    dim = 4;
}

Matrix::~Matrix() {}

void Matrix::multiply(double p[]) {
    double result[4] = {0,0,0};
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 3; j++) {
            result[i] += values[i][j] * p[j];
        }
    }
    for(int i = 0; i < 3; i++)
        p[i] = result[i];
}