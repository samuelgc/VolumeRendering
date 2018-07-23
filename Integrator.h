#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <vector>
#include <cmath>

#include "Volume.h"

using namespace std;

class Integrator {

public:

    Integrator();

    ~Integrator();

    /**
     * Intersects a ray with a given bounding box
     * 
     * @param double orig[] - ray origin
     * @param double d[] - ray direction
     * @param Volume* box - bounding box to intersect against
     * @param double *t_far - exiting intersecinot with bbox
     * @return the nearest point of intersection on bbox surface
     */
    double intersect(double orig[], double d[], Volume* box, double *t_far);

    /**
     * Performs a raymarch for the given ray from the orig in direction d
     * Stores the resulting color in result
     * 
     * @param double orig[] - ray origin
     * @param double d[] - ray direction
     * @param vector<Volume*> objs - objects to sample from
     * @param double result[] - resulting color value from integration
     */
    void integrate(double orig[], double d[], vector<Volume*> objs, double result[]);

    /**
     * Computes the radiance at the given position with the incoming
     * ray specified
     * 
     * @param double pos[] - position to calculate radiance for
     * @param double dir[] - direction of incoming ray
     * @param Volume* obj - the volume to sample radiance from
     * @param double rgb[] - resulting color value
     */ 
    void radiance(double pos[], double dir[], Volume* obj, double rgb[]);

};


#endif //INTEGRATOR_H
