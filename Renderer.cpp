#include "Renderer.h"

Renderer::Renderer() {}

Renderer::~Renderer() {}

void Renderer::loadScene(string inFile) {
    scene = new Scene(inFile);
}

void Renderer::render() {

}

void Renderer::write(string outFile) {

}