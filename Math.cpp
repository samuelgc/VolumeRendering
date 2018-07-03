#include "Math.h"


double dist(double a[], double b[]){
    double sum = 0;
    for (int i = 0; i < 3; i++)
        sum += pow(b[i] - a[i], 2.0);
    return sqrt(sum);
}

Matrix* mat_mult(Matrix* a, Matrix* b) {
    Matrix* result = new Matrix();
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            double value = 0;
            for (int k = 0; k < 4; k++)
                value += a->get(i, k) * b->get(k, j);
            result->set(i, j, value);
        }
    }
    return result;
}

int convert(double value) {
    value *= 255;
    if(value > 255)
        return 255;
    if(value < 0)
        return 0;
    return int(value);
}