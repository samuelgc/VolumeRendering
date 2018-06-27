#ifndef CAMERA_H
#define CAMERA_H


class Camera {

public:

    Camera();

    Camera(double lookFrom[], double lookAt[], int x, int y);

    ~Camera();

    /**
     * returns the look from point of the Camera
     */
    double* getOrigin();

    /**
     * returns the look at point of the Camera
     */
    double* getOrient();

    /**
     * returns the resolution of the screen as an array { x, y }
     */
    int* getRes();

private:

    double origin[3];
    double orient[3];
    int res[2];

};


#endif CAMERA_H
