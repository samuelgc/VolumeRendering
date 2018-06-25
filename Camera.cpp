#include "Camera.h"

Camera::Camera() {
    for (int i = 0; i < 3; i++) {
        origin[i] = 0;
        orient[i] = 0;
        if (i < 2)
            res[i] = 100;
    }
    origin[2] = 10;
}

Camera::Camera(double lookFrom[], double lookAt[], int x, int y) {
    for (int i = 0; i < 3; i++) {
        origin[i] = lookFrom[i];
        orient[i] = lookAt[i];
    }
    res[0] = x;
    res[1] = y;
}

Camera::~Camera() {}

double* Camera::getOrigin() {
    return origin;
}

double* Camera::getOrient() {
    return orient;
}

int* Camera::getRes() {
    return res;
}