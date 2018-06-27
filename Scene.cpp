#include "Scene.h"

#include <iostream>

Scene::Scene() {
    cam = new Camera();
}

Scene::Scene(string filename) {
    parse(filename);
}

Scene::~Scene() {}

void Scene::parse(string filename) {
    ifstream file;
    file.open(filename);
    if (!file.is_open()) {
        cout << "Failed to load scene.\n";
        return;
    }

    // Parse out elements of the scene

    file.close();
}

int Scene::width() {
    return cam->getRes()[0];
}

int Scene::height() {
    return cam->getRes()[1];
}

void Scene::setTransform() {

}