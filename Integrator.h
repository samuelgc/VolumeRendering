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
     * @return the nearest point of intersection on bbox surface
     */
    double intersect(double orig[], double d[], Volume* box);

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

};


#endif //INTEGRATOR_H
