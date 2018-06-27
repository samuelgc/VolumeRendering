#include "Renderer.h"

Renderer::Renderer() {}

Renderer::~Renderer() {}

void Renderer::loadScene(string inFile) {
    scene = new Scene(inFile);
}

void Renderer::render() {
    // Assume a camera field of view of 30 degrees

}

void Renderer::write(string outFile) {
    ofstream file;
    string filename = "..\\results\\" + outFile;
    file.open(filename);
    file << "P3\n";
    file << scene->width() << " " << scene->height() << "\n";
    file << "255\n";
    for (int i = 0; i < pixels.size(); i++)
        file << (int)pixels.at(i) << " ";
    file.close();
}