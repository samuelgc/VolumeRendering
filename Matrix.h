#ifndef MATRIX_H
#define MATRIX_H

/**
 * Class of 4x4 Matrices
 */
class Matrix {


public:

    Matrix();

    ~Matrix();

    /**
     * Set all values to 0 except {3, 3} which is set to 1
     * This is for use in homogeneous coordinate transformations
     */
    void reset();

    /**
     * Sets the value at {i, j} to the value given
     *
     * @param i
     * @param j
     * @param value
     */
    void set(int i, int j, double value);

    /**
     * Returns the value at location {i, j}
     *
     * @param i
     * @param j
     * @return
     */
    double get(int i, int j);

    /**
     * Multiplies the given vector p by this Matrix
     *
     * @param p
     */
    void multiply(double p[]);

private:

    double values[4][4];

};


#endif //MATRIX_H
