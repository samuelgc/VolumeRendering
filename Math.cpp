#include "Math.h"


double dist(double a[], double b[]){
    double sum = 0;
    for (int i = 0; i < 3; i++)
        sum += pow(b[i] - a[i], 2.0);
    return sqrt(sum);
}

double len(double a[]) {
    double length = 0;
    for(int i = 0; i < 3; i++)
        length += a[i] * a[i];
    return sqrt(length);
}

void sum(double a[], double b[]) {
    for(int i = 0; i < 3; i++)
        a[i] += b[i];
}

void subtract(double a[], double b[]) {
    for(int i = 0; i < 3; i++)
        b[i] = b[i] - a[i];
}

void scale(double a[], double s) {
    for(int i = 0; i < 3; i++) 
        a[i] *= s;
}

void normalize(double a[]) {
    double div = len(a);
    for(int i = 0; i < 3; i++)
        a[i] /= div;
}

void move(double o[], double d[], double t, double res[]) {
    for(int i = 0; i < 3; i++)
        res[i] = o[i] + d[i] * t;
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