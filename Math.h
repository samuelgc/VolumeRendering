#ifndef MATH_H
#define MATH_H


#include <cmath>
#include "Matrix.h"

const double PI = 3.14159265359;

/// Vector Math ///

/**
 * Calculates the distance between vectors a and b
 *
 * @param a - vector
 * @param b - vector
 * @return
 */
double dist(double a[], double b[]);

/**
 * Subtracts a from b (i.e. b - a)
 * 
 * @param a - vector
 * @param b - vector
 */
void subtract(double a[], double b[]);

/// Matrix Math ///

/**
 * Multiplies two matrices and returns the resulting Matrix
 *
 * @param a
 * @param b
 * @return
 */
Matrix* mat_mult(Matrix* a, Matrix* b);

/// Color Math ///

/**
 * Converts a float 0-1 color to int 0-255 color
 *
 * @param value
 * @return
 */
int convert(double value);

#endif //MATH_H
