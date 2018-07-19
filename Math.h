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

/**
 * Scales vector a by value s
 * 
 * @param a - vector
 * @param s - scalar value
 */
void scale(double a[], double s);

/**
 * Moves a point along parametrized ray to location at t
 *
 * @param o - origin
 * @param d - direction of ray
 * @param t - time
 * @param res - resulting location
 */
void move(double o[], double d[], double t, double res[]);

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
