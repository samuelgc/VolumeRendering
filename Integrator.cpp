#include "Integrator.h"

Integrator::Integrator() {}

Integrator::~Integrator(){}

double Integrator::intersect(double orig[], double d[], Volume* box) {
    double t_near = -INFINITY;
    double t_far = INFINITY;
    double t1, t2;
    for(int i = 0; i < 3; i++) {
        if(d[i] == 0) {
            if(orig[i] < box->getMin()[i] || orig[i] > box->getMax()[i])
                return -1;
        }
        t1 = (box->getMin()[i] - orig[i]) / d[i];
        t2 = (box->getMax()[i] - orig[i]) / d[i];
        if(t1 > t2) {
            double temp = t1;
            t1 = t2;
            t2 = temp;
        }
        if(t1 > t_near)
            t_near = t1;
        if(t2 < t_far)
            t_far = t2;
        if(t_near > t_far)
            return -1;
        if(t_far < 0)
            return -1;
    }
    return t_near;
}

void Integrator::integrate(double orig[], double d[], vector<Volume*> objs, double result[]) {
    // TODO: Here's the part where we shoot the ray through and get the color back
    vector<Volume*> hit;
    for(Volume* vol : objs) {
        if(intersect(orig, d, vol) != -1)
            hit.push_back(vol);
    }
    if(!hit.empty()) {
        // Not sure how to handle hitting multiple volumes actually... need to think about it
    }
}

double* Integrator::radiance(double pos[], double dir[], Volume* v) {
    Material* m = v->getMat();

    double density = m->dense_scale() * v->sample(pos, 0); // Sample density?
    double* smokecolor = m->dense_color();
    for(int i = 0; i < 3; i++)
        smokecolor[i] *= m->dense_intense();
    
    double emit = m->temp_intense() * v->sample(pos, 1); // Sample heat?
    double temp = m->fire_intense() * v->sample(pos, 2); // Sample temperature?
    temp *= m->kelvin_temp();

    // Perform blackbody mapping from temperature to RGB
}
