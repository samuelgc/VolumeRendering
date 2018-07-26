#ifndef RENDERER_H
#define RENDERER_H


#include <fstream>
#include <string>
#include <vector>

#include "Scene.h"
#include "Integrator.h"

using namespace std;

class Renderer {

public:

    Renderer();

    ~Renderer();

    /**
     * Creates an empty default scene
     */
    void newScene();

    /**
     * Loads a scene descritption from a file
     */
    void loadScene(string inFile);

    /**
     * Renders the scene to the pixel buffer 
     * 
     * @param samples - ray samples
     */
    void render(double samples);

    /**
     * Writes the pixel buffer out to an img file
     */
    void write(string outFile);

private:

    /**
     * Converts float rgb values to int and adds to pixel buffer
     *
     * @param rgb
     */
    void addPixel(double rgb[]);

    Scene* scene;
    vector<int> pixels;

    Integrator* march;

};


#endif //RENDERER_H
