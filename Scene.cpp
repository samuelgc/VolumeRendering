#include "Scene.h"

#include <iostream>
#include <cmath>

#include "Math.h"
#include "dataGetter/loadGeo2.0.h"

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
    string waste;
    file >> waste;              /// READ CAMERA
    double lf[3] = {0,0,0};
    double la[3] = {0,0,0};
    int resx, resy;
    for(int i = 0; i < 3; i++)
        file >> lf[i];
    for(int i = 0; i < 3; i++)
        file >> la[i];
    file >> resx;
    file >> resy;
    cam = new Camera(lf, la, resx, resy);

    // TODO: potentially add lights
    
    // TODO: Read in Volumes
    string fire_file;
    file >> fire_file;
    send_vol_data svd = getAllData(fire_file);
    Volume* v = new Volume();
    v->loadFireData(svd);
    addVolume(v);

    file.close();
}

int Scene::width() {
    return cam->getRes()[0];
}

int Scene::height() {
    return cam->getRes()[1];
}

double* Scene::origin() {
    return cam->getOrigin();
}

void Scene::setCorner() {
    // Set size of screen in world units
    double len = dist(cam->getOrigin(), cam->getOrient());
    len *= tan(.4);// Assume a camera field of view of 45 degrees roughly...
    corner[0] = -1 * len;
    corner[1] = len * (cam->getRes()[1] / cam->getRes()[0]);
    inc = (len * 2.0) / cam->getRes()[0];
}

double* Scene::getCorner() {
    return corner;
}

double Scene::getInc() {
    return inc;
}

void Scene::setTransform() {
    transfo = new Matrix();
    for(int i = 0; i < 4; i++)
        transfo->set(i, i, 1);
    /*
    Matrix* mult = new Matrix();
    for(int i = 0; i < 4; i++)
        mult->set(i, i, 1);
    double point[3] = {0,0,0};
    for(int i = 0; i < 3; i++)
        point[i] = cam->getOrigin()[i] - cam->getOrient()[i];
    
    // Rotate around Y --- THIS MIGHT BE WRONG BECAUSE OF DIRECTION OF Z
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
    */
    transfo->set(0, 3, cam->getOrient()[0]);
    transfo->set(1, 3, cam->getOrient()[1]);
    transfo->set(2, 3, cam->getOrient()[2]);

    setCorner();
}

void Scene::transform(double p[]) {
    transfo->multiply(p);
}

void Scene::addVolume(Volume* v) {
    volumes.push_back(v);
}

vector<Volume*> Scene::getVolumes() {
    return volumes;
}
