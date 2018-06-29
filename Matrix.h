#ifndef MATRIX_H
#define MATRIX_H


class Matrix {


public:

    Matrix();

    ~Matrix();

    void multiply(double p[]);

private:

    double values[4][4];
    int dim;

};


#endif //MATRIX_H
