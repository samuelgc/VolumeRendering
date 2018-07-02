#ifndef MATH_H
#define MATH_H


#include <cmath>
#include "Matrix.h"

const double PI = 3.14159265359;

/**
 * Calculates the distance between vectors a and b
 *
 * @param a - vector
 * @param b - vector
 * @return
 */
double dist(double a[], double b[]);

/**
 * Multiplies two matrices and returns the resulting Matrix
 *
 * @param a
 * @param b
 * @return
 */
Matrix* mat_mult(Matrix* a, Matrix* b);

#endif //MATH_H
