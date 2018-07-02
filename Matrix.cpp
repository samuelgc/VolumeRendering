#include "Matrix.h"


Matrix::Matrix() {
    reset();
}

Matrix::~Matrix() {}

void Matrix::reset() {
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            values[i][j] = 0;
        }
    }
}

void Matrix::set(int i , int j, double value) {
    values[i][j] = value;
}

double Matrix::get(int i, int j) {
    return values[i][j];
}

void Matrix::multiply(double p[]) {
    double result[3] = {0,0,0};
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            result[i] += values[i][j] * p[j];
        }
        result[i] += values[i][3];
    }
    for(int i = 0; i < 3; i++)
        p[i] = result[i];
}