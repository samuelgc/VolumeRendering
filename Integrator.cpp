#include "Integrator.h"

#include <stdlib.h>

#include "Blackbody.h"
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

/**
 * Samuel version...
 *
void Integrator::integrate(double orig[], double d[], vector<Volume*> objs, double result[]) {
    Volume* vol = objs.at(0);
    double eor = 0;
    double t = intersect(orig, d, vol, &eor);
    if(t < 0) {
        return;
    }

    // Move intersections to edges of the volume
    double pos[3] = {0, 0, 0};
    move(orig, d, eor, pos);
    while(vol->sample(pos, 0) <= 0.001) {
        eor -= vol->getSize();
        move(orig, d, eor, pos);
        if(eor < t)
            return;
    }
    do {
        move(orig, d, t, pos);
        t += vol->getSize();
        if(t > eor)
            return;
    } while(vol->sample(pos, 0) <= 0.001);
    t -= vol->getSize();

    // Perform path trace
    double wig = 0;
    double w[3] = {d[0], d[1], d[2]};
    eor -= t;
    while(true) {
        wig = (double)rand() / RAND_MAX;
        if(wig > eor)
            return;
        move(pos, w, wig, pos);
        if(wig < .333) {
            double rgb[3] = {1,1,1};
            radiance(pos, w, vol, rgb);
            sum(result, rgb);
            return;
        } else if (wig < 1 - .333) {
            eor = 0;
            for(int i = 0; i < 3; i++)
                w[i] = (double)rand() / RAND_MAX -.5;
            normalize(w);
            double temp[3] = {pos[0], pos[1], pos[2]};
            while(vol->sample(temp, 0) > 0.001) {
                eor += vol->getSize();
                move(pos, w, eor, temp);
            }
            eor -= vol->getSize();
        } else
            eor -= wig;
    }
}

/**
 * Brian version...
 */
void Integrator::integrate(double orig[], double d[], vector<Volume*> objs, double result[]) {
    Volume* vol = objs.at(0);
    double eor = 0;
    double t = intersect(orig, d, vol, &eor);
    if(t < 0) {
        return;
    }

    // Move intersections to edges of the volume
    double pos[3] = {0, 0, 0};
    move(orig, d, eor, pos);
    while(vol->sample(pos, 0) <= 0.0001) 
    {   
        eor -= vol->getSize();
        move(orig, d, eor, pos);
        if(eor < t)
            return;
    }
    move(orig,d,t,pos);
    while(vol->sample(pos, 0) <= 0.0001)
    {
        if(t > eor){
            return;
        }
        t+= vol->getSize();
        move(orig,d,t,pos);
    }
    // t -= vol->getSize();
    double absor = 1.0;
    double scat = 3.0;
    double exti = absor + scat;
    double nc = 5;
    double maj = exti + nc; 
    double ray_length = eor - t;
    // Perform path trace
    // double wig = 0;
    double w[3] = {d[0], d[1], d[2]};
    // eor -= t;
    // this might not be needed
    move(orig, d, t, pos);
    while(true)
    {
        double r1 = (double)rand() / RAND_MAX;
        double r2 = (double)rand() / RAND_MAX;
        double dis = abs(log(1-r1)/maj);
        // cout << dis << endl;
        if(dis > ray_length) // the next move takes it out of the volumes(fire) boundary has been hit
        {
            // double rgb[3] = {0,255,255}; // for debugging
            // sum(result,rgb);
            return;
        }
        //move particle to next position
        move(pos,w,dis,pos);
        if(r2 < absor/maj)
        {
            // cout << vol->sample(pos,0) << endl;
            double rgb[3] = {1,1,1};//{50/255.0,255/255.0,50/255.0};
            radiance(pos, w, vol, rgb);
            sum(result, rgb);
            return;
        }
        // else if(r2 < 1 - nc/maj)
        // {
        //     // double rgb[3] = {0,0,0};
        //     // double rgb[3] = {50/255.0,255.0/255.0,50/255.0};//{50,250,50};
        //     // double rgb[3] = {0,0,0};
        //     for(int i = 0 ; i < 3 ; i++)
        //         w[i] = (double)rand() / RAND_MAX;
        //     normalize(w);
        //     double temp[3] = {pos[0], pos[1], pos[2]};
        //     ray_length = 0;
        //     int count = 0;
        //     while(vol->sample(temp,0) < .0001 )
        //     {
        //         if(count > 100)
        //             cout << "im stuck" << endl;
        //         ray_length += vol->getSize();
        //         move(pos,w,ray_length,temp);
        //     }
        //     // sum(result, rgb);
        //     // return;
        // }
        else
        {
            ray_length -= dis;
        }
    }
}

void Integrator::radiance(double pos[], double dir[], Volume* v, double rgb[]) {
    Material* m = v->getMat();
   
    double density = m->dense_scale() * v->sample(pos, 0); // Sample density?
    double emit = m->temp_intense() * v->sample(pos, 5); // Sample heat?
    
    if(emit == 0) {
        for(int i = 0; i < 3; i++)
            // rgb[i] = 0;
            rgb[i] = m->dense_color()[i] * m->dense_intense();
        scale(rgb, density);
        return;
    }
    double temp = m->fire_intense() * v->sample(pos, 4); // Sample temperature?
    temp *= m->kelvin_temp();

    
    // Perform blackbody mapping from temperature to XYZ
    double chr0[3] = {0,0,0};
    double chr1[3] = {0,0,0};
    double chr2[3] = {0,0,0};
    double adapt_val = m->adaption();
    double burn_val = m->burn();
    double val_x = blackbody(temp, adapt_val, burn_val, chr0);
    double val_y = blackbody(temp, adapt_val, burn_val, chr1);
    double val_z = blackbody(temp, adapt_val, burn_val, chr2);

    double xyz[3] = {0,0,0};
    xyz[0] = chr0[0] * val_x;
    xyz[1] = chr1[1] * val_y;
    xyz[2] = chr2[2] * val_z;

    //cout << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << " => ";

    // Convert XYZ to sRGB colorspace
    rgb[0] = 3.2404542 * xyz[0] - 1.5371385 * xyz[1] - 0.4985314 * xyz[2];
    rgb[1] = -0.969266 * xyz[0] + 1.8760108 * xyz[1] + 0.041556 * xyz[2];
    rgb[2] = 0.0556434 * xyz[0] - 0.2040259 * xyz[1] + 1.0572252 * xyz[2];
    for(int i = 0; i < 3; i++) {
        if(rgb[i] <= 0.0031308)
            rgb[i] *= 12.92;
        else
            rgb[i] = pow((rgb[i] * 1.055), (1.0 / 2.4)) - 0.055;
    }

    //cout << rgb[0] << ", " << rgb[1] << ", " << rgb[2] << "\n";
}
