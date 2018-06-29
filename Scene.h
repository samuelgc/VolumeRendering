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
     * Set screen to world transform
     */
    void setTransform();

    /**
     * Transforms a point from screen space to world space
     *
     * @param p
     */
    void transform(double p[]);

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

};


#endif //SCENE_H
