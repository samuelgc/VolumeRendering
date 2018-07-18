#include "Material.h"


Material::Material() { 
    d_scale = 1;
    d_int = 1;
    for(int i = 0; i < 3; i++)
        d_rgb[i] = 0.2;

    f_int = 2;
    t_int = .2;
    kelvin = 5000;
    a = 0.15;
    b = 0;     
}

Material::~Material() {}

void Material::setAll(vector<double> values) {
    d_scale = values[0];
    d_int = values[1];
    for(int i = 0; i < 3; i++)
        d_rgb[i] = values[i+2];

    f_int = values[5];
    t_int = values[6];
    kelvin = values[7];
    a = values[8];
    b = values[9];
}

double Material::dense_scale() {
    return d_scale;
}

double Material::dense_intense() {
    return d_int;
}

double* Material::dense_color() {
    return d_rgb;
}

double Material::fire_intense() {
    return f_int;
}

double Material::temp_intense() {
    return t_int;
}

double Material::kelvin_temp() {
    return kelvin;
}

double Material::adaption() {
    return a;
}

double Material::burn() {
   return b;
}