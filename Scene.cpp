#include "Scene.h"
#include "Math.h"

#include <iostream>
#include <cmath>

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
    transfo = new Matrix();
    transfo->set(3, 3, 1);
    Matrix* mult = new Matrix();
    mult->set(3, 3, 1);
    double point[3] = {0,0,0};
    for(int i = 0; i < 3; i++)
        point[i] = cam->getOrigin()[i] - cam->getOrient()[i];

    // Rotate around Y
    double theta = -1.0 * (atan2(point[2], point[0]) * (180.0 / PI) + 90.0);
    double theta_x = -1.0 * theta;
    double sin_t = sin(theta);
    double cos_t = cos(theta);
    mult->set(0, 0, cos_t);
    mult->set(0, 2, sin_t);
    mult->set(2, 0, -sin_t);
    mult->set(2, 2, cos_t);
    mult->multiply(point);

    // Rotate around X
    theta = atan2(point[1], point[2]) * (180.0 / PI) + 90.0;
    sin_t = sin(theta);
    cos_t = cos(theta);
    transfo->set(1, 1, cos_t);
    transfo->set(1, 2, -sin_t);
    transfo->set(2, 1, sin_t);
    transfo->set(2, 2, cos_t);

    mult->reset();
    sin_t = sin(theta_x);
    cos_t = cos(theta_x);
    mult->set(0, 0, cos_t);
    mult->set(0, 2, sin_t);
    mult->set(2, 0, -sin_t);
    mult->set(2, 2, cos_t);
    transfo = mat_mult(mult, transfo);

    // Translate
    mult->reset();
    for(int i = 0; i < 3; i++)
        mult->set(i, i, 1);
    mult->set(0, 3, cam->getOrient()[0]);
    mult->set(1, 3, cam->getOrient()[1]);
    mult->set(2, 3, cam->getOrient()[2]);
    transfo = mat_mult(mult, transfo);
}

void Scene::transform(double p[]) {
    transfo->multiply(p);
}
