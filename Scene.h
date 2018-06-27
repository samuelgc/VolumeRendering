#include <fstream>
#include <string>
#include <vector>

#include "Camera.h"
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

private:

    /**
     * Reads a file and populates Scene parameters
     * * Sets up Camera
     * * Adds Volumes to scene
     */
    void parse(string filename);

    Camera* cam;
    vector<Volume*> volumes;

};