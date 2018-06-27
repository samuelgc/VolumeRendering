#include "Math.h"


double dist(double a[], double b[]){
    double sum = 0;
    for (int i = 0; i < 3; i++)
        sum += pow(b[i] - a[i], 2.0);
    return sqrt(sum);
}
