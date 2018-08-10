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
    if(t < 0) {
        return;
    }

    // Move intersections to edges of the volume
    double pos[3] = {0, 0, 0};
    move(orig, d, eor, pos);
    while(vol->sample(pos, 5) <= 0.0001) 
    {   
        // if(debug)
        //     cout << pos[0] << "," << pos[1] << "," << pos[2] << " end of array " << eor << endl;
        eor -= vol->getSize();
        move(orig, d, eor, pos);
        if(eor < t)
            return;
    }
    move(orig,d,t,pos);
    while(vol->sample(pos, 5) <= 0.0001)
    {
        if(t > eor){
            return;
        }
        t+= vol->getSize();
        move(orig,d,t,pos);
    }
    // do {
        
    //     move(orig, d, t, pos);
    //     t += vol->getSize();
    //     // if(debug)
    //     //     cout << pos[0] << "," << pos[1] << "," << pos[2] << " end of array " << eor << endl;
    //     if(t > eor)
    //         return;
    // } while(vol->sample(pos, 0) <= 0.0001);
    // t -= vol->getSize();
    // if(debug)
        // cout << t << " " << eor << endl;
    // double rgb[3] = {0,0,0};
    // sum(result, rgb);
    //return;

    double absor = 1.0;
    double scat = 1.0;
    double exti = absor + scat;
    double nc = 2;
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
        if(true)//r2 < absor/maj)
        {
            // cout << vol->sample(pos,0) << endl;
            double rgb[3] = {1,1,1};//{50/255.0,255/255.0,50/255.0};
            radiance(pos, w, vol, rgb);
            sum(result, rgb);
            return;
        }
        else if(r2 < 1 - nc/maj)
        {
            double rgb[3] = {0,0,0};
            // double rgb[3] = {0,255.0/255,255.0/255};//{50,250,50};
            sum(result, rgb);
            return;
        }
        else
        {
            double rgb[3] = {0,0,0};
            // double rgb[3] = {255.0/255,255.0/255,0};//{50,250,50};
            sum(result, rgb);
            return;
        }
        
        // return;

    }
    // while(true) {
    //     wig = (double)rand() / RAND_MAX;
    //     if(wig > eor)
    //     {
            
    //         return;
    //     }
    //     move(pos, w, wig , pos); // wig
    //     if(true){//wig < .6) { // was .6
    //         double rgb[3] = {1,1,1};
    //         radiance(pos, w, vol, rgb);
    //         if(debug)
    //             cout << rgb[0] << " " << rgb[1] << " " << rgb[2] << endl; 
    //         sum(result, rgb);
    //         return;
    //     } else if (wig < 1 - .2) {
    //         eor = 0;
    //         for(int i = 0; i < 3; i++)
    //             w[i] = (double)rand() / RAND_MAX -.5;
    //         normalize(w);
    //         double temp[3] = {pos[0], pos[1], pos[2]};
    //         while(vol->sample(temp, 0) > 0.001) {
    //             eor += vol->getSize();
    //             move(pos, w, eor, temp);
    //         }
    //         eor -= vol->getSize();
    //     } else
    //     {
    //         eor -= wig;
    //     }
    // }
}

void Integrator::radiance(double pos[], double dir[], Volume* v, double rgb[]) {
    Material* m = v->getMat();
    // while(v->sample(pos,5) == 0)
    // {
    //     move(pos,dir,.05,pos);
    // }
    double density = m->dense_scale() * v->sample(pos, 0); // Sample density?
    double emit = m->temp_intense() * v->sample(pos, 5); // Sample heat?
    
    if(emit == 0) {
        for(int i = 0; i < 3; i++)
            rgb[i] = 0;
            // rgb[i] = m->dense_color()[i] * m->dense_intense();
        // scale(rgb, density*5);
        return;
    }
    double temp = m->fire_intense() * v->sample(pos, 4); // Sample temperature?
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

    scale(rgb, 0.00392156863); // divide by 255
    if(rgb[0] < 0 || rgb[0] > 255 || rgb[1] < 0 || rgb[1] > 255 || rgb[2] < 0 || rgb[2] > 255)
        cout << "Out of Range\n";
    // scale(rgb, density*5); // Do this?
}
