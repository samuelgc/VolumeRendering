#include "Integrator.h"

#include <stdlib.h>

#include "Math.h"

Integrator::Integrator() {}

Integrator::~Integrator(){}

double Integrator::intersect(double orig[], double d[], Volume* box, double *t_far) {
    double t_near = -INFINITY;
    *t_far = INFINITY;
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
        if(t2 < *t_far)
            *t_far = t2;
        if(t_near > *t_far)
            return -1;
        if(*t_far < 0)
            return -1;
    }
    return t_near;
}

void Integrator::integrate(double orig[], double d[], vector<Volume*> objs, double result[]) {
    Volume* vol = objs.at(0);
    double eor = 0;
    double t = intersect(orig, d, vol, &eor);
    if(t == -1)
        return;
    double pos[3] = {orig[0], orig[1], orig[2]};
    double density = 0;
    while(density <= 0) {
        t += vol->getSize();
        move(orig, d, t, pos);
        density = vol->sample(pos, 0);
        if(t > eor)
            return;
    }

    double wig = 0;
    double w[3] = {d[0], d[1], d[2]};
    while(true) {
        wig = (double)rand() / RAND_MAX;
        // t = -1.0 * log(1.0 - wig);
        // if(t > end_of_ray)
        //      return;
        // subtract(pos, scaled(w, t));
        // if(wig < absorption_probablity_at_current_location) {
        //     sum(result, radiance(pos, w, hit.at(0)));
        //     return;
        // } else if (wig < 1 - scattering_probability_at_current_location)
        //      update w based on phase function (isotropic)
        //      update the end_of_ray
        // } else
        //      end_of_ray = end_of_ray - t
    }
}

double* Integrator::radiance(double pos[], double dir[], Volume* v) {
    Material* m = v->getMat();
    double rgb[3] = {0,0,0};

    double density = m->dense_scale() * v->sample(pos, 0); // Sample density?
    double emit = m->temp_intense() * v->sample(pos, 1); // Sample heat?
    if(emit == 0) {
        for(int i = 0; i < 3; i++)
            rgb[i] = m->dense_color()[i] * m->dense_intense();
        scale(rgb, density);
        return rgb;
    }

    double temp = m->fire_intense() * v->sample(pos, 2); // Sample temperature?
    temp *= m->kelvin_temp();

    // Perform blackbody mapping from temperature to RGB
    temp /= 100;
    if(temp <= 66) {
        rgb[0] = 255;

        rgb[1] = temp;
        rgb[1] = 99.4708025861 * log(rgb[1]) - 161.1195681661;

        if(temp <= 19)
            rgb[2] = 0;
        else {
            rgb[2] = temp - 10;
            rgb[2] = 138.5177312231 * log(rgb[2]) - 305.0447927307;
            if(rgb[2] < 0)
                rgb[2] = 0;
            if(rgb[2] > 255)
                rgb[2] = 255;
        }
    } else {
        rgb[0] = temp - 60;
        rgb[0] = 329.698727446 * pow(rgb[0], -0.1332047592);
        if(rgb[0] < 0)
            rgb[0] = 0;
        if(rgb[0] > 255)
            rgb[0] = 255;
        
        rgb[1] = temp - 60;
        rgb[1] = 288.1221695283 * pow(rgb[1], -0.0755148492);
        rgb[2] = 255;
    }
    if(rgb[1] < 0)
        rgb[1] = 0;
    if(rgb[1] > 255)
        rgb[1] = 255;

    scale(rgb, 0.00392156863);
    scale(rgb, density); // Do this?
    return rgb;
}
