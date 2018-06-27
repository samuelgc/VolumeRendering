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

    int width();

    int height();

private:

    void parse(string filename);

    Camera* cam;
    vector<Volume*> volumes;

};