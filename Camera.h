

using namespace std;

class Camera {

public:

    Camera();

    Camera(double lookFrom[], double lookAt[], int x, int y);

    ~Camera();

    double* getOrigin();

    double* getOrient();

    int* getRes();

private:

    double origin[3];
    double orient[3];
    int res[2];

};