#ifndef SCENE_H
#define SCENE_H


#include <fstream>
#include <string>
#include <vector>

#include "Camera.h"
#include "Matrix.h"
#include "Volume.h"

using namespace std;

class Scene {

public:

    Scene();

    /**
     * Loads a scene from file description
     */
    Scene(string filename);

    ~Scene();

    /**
     * returns the X frame resolution
     */
    int width();

    /**
     * returns the Y frame resolution
     */
    int height();

    /**
     * returns the location of the camera
     */
    double* origin();

    /**
     * Sets world coordinate of screen top left corner
     */
    void setCorner();

    /**
     * return top left hand corner in world space
     */
    double* getCorner();

    /**
     * returns the size of a pixel
     */
    double getInc();

    /**
     * Set screen to world transform
     */
    void setTransform();

    /**
     * Transforms a point from screen space to world space
     *
     * @param p
     */
    void transform(double p[]);

    /**
     * Adds a volume to the scene
     * 
     * @param v - volume to add
     */
    void addVolume(Volume* v);

    /**
     * Returns a list of all the volumes in the scene
     */
    vector<Volume*> getVolumes();

private:

    /**
     * Reads a file and populates Scene parameters
     * * Sets up Camera
     * * Adds Volumes to scene
     */
    void parse(string filename);

    Camera* cam;
    vector<Volume*> volumes;

    Matrix* transfo;
    double corner[2];
    double inc;
};


#endif //SCENE_H
