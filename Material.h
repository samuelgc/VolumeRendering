#ifndef MATERIAL_H
#define MATERIAL_H

#include <vector>

using namespace std;

class Material {

public:

    /**
     * Default Material constructor
     */
    Material();

    ~Material();

    /**
     * Sets all the private variables in order of listing
     * 
     * @param vector<double> values - ordered values to set
     */
    void setAll(vector<double> values);

    /**
     * Returns the Density Scale
     */
    double dense_scale();

    /**
     * Returns the Smoke Brightness
     */
    double dense_intense();

    /**
     * Returns the Smoke Color
     */
    double* dense_color();

    /**
     * Returns the Fire Intensity
     * which acts as a scalar on the heat field
     */
    double fire_intense();

    /**
     * Returns the Temperature Intensity
     * which acts as a scalar on the temperature field
     */
    double temp_intense();

    /**
     * Returns the Color Temperature (in Kelvin)
     * which is used to map temperature to color values
     */
    double kelvin_temp();

    /**
     * Returns the Adaption Value
     */
    double adaption();

    /**
     * Returns the Burn value
     */
    double burn();


private:

    // Smoke parameters
    double d_scale;     // Density Scale
    double d_int;       // Smoke Brightness
    double d_rgb[3];    // Smoke Color
    
    // Fire parameters
    double f_int;       // Intensity Scale
    double t_int;       // Temperature Scale
    double kelvin;      // Color Temperature
    double a;           // Adaption
    double b;           // Burn

};


#endif // MATERIAL_H 