#include "Renderer.h"

#include "Math.h"

Renderer::Renderer() {
    march = new Integrator;
}

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
    double point[3] = {0, 0, 0};
    double rgb[3] = {0,0,0};
    double old_y = scene->getCorner()[1];
    for(int x = 0; x < scene->width(); x++) {
        double old_x = scene->getCorner()[0];
        for(int y = 0; y < scene->height(); y++) {
            point[2] = 0;
            point[1] = old_y;
            point[0] = old_x;
            scene->transform(point);
            subtract(scene->origin(), point);
            for(int i = 0; i < 64; i++)
                march->integrate(scene->origin(), point, scene->getVolumes(), rgb);
            scale(rgb, 0.015625);
            addPixel(rgb);
            old_x += i;
        }
        old_y -= i;
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