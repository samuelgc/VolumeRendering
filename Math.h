#ifndef MATH_H
#define MATH_H


#include <cmath>
#include "Matrix.h"

const double PI = 3.14159265359;

/**
 * Returns a if a > b else returns b
 */
double max(double a, double b);

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
 * Return the length of vector a
 *
 * @param a - vector
 * @return
 */
double len(double a[]);

/**
 * Adds the vector b to a
 *
 * @param a - vector
 * @param b - vector
 */
void sum(double a[], double b[]);

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
 * Divides vector a by value div
 * 
 * @param a - vector
 * @param div - dividend
 */
void divide(double a[], double div);

/**
 * Normalizes the vector given
 *
 * @param a - vector
 */
void normalize(double a[]);

/**
 * Copies values in a to b
 * 
 * @param a - vector to copy
 * @param b - destination vector
 */
void copy(double a[], double b[]);

/**
 * Sets all values of a to 0
 */
void reset(double a[]);

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
