#include "Renderer.h"

#include "Math.h"

Renderer::Renderer() {}

Renderer::~Renderer() {}

void Renderer::newScene() {
    scene = new Scene();
}

void Renderer::loadScene(string inFile) {
    scene = new Scene(inFile);
}

void Renderer::render() {
    // Assume a camera field of view of 30 degrees
    scene->setTransform();
    double i = scene->getInc();
    double point[3] = {0,0,0};
    double rgb[3] = {0,0,0};
    for(int x = 0; x < scene->width(); x++) {
        point[1] = 0;
        for(int y = 0; y < scene->height(); y++) {
            point[2] = 0;
            point[0] += i;
            point[1] += i;
            scene->transform(point);
            // TODO: Here's the part where we shoot the ray through and get the color back
            addPixel(rgb);
        }
    }
}

void Renderer::write(string outFile) {
    ofstream file;
    string filename = "output\\" + outFile;
    file.open(filename);
    file << "P3\n";
    file << scene->width() << " " << scene->height() << "\n";
    file << "255\n";
    for(int pixel : pixels)
        file << pixel << " ";
    file.close();
}

void Renderer::addPixel(double rgb[]) {
    for(int i = 0; i < 3; i++)
        pixels.push_back(convert(rgb[i]));
}