#include <fstream>
#include <string>
#include <vector>

#include "Scene.h"

using namespace std;

class Renderer {

public:

    Renderer();

    ~Renderer();

    /**
     * Loads a scene descritption from a file
     */
    void loadScene(string inFile);

    /**
     * Renders the scene to the pixel buffer 
     */
    void render();

    /**
     * Writes the pixel buffer out to an img file
     */
    void write(string outFile);

private:

    Scene* scene;
    vector<int> pixels;

};