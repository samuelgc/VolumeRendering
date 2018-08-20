#ifndef BLACKBODY_H
#define BLACKBODY_H

#define C_PI                3.14159265358979323846  // pi
#define PLANCK_S            5.67037321e-8           // Stefan-Boltzmann constant (W/m^2)
#define PLANCK_LPA_D10      6.94559                 // = planckbandlogpoweravg(0,1000)
#define PLANCK_TMNORM       568.0                   // unit^-2


/** 
 * The total spectral radiant power (W/m^2) for a surface (i.e: total 
 * hemispherical emittance incl. Lambert's cos term) at temperature T  
 * for lambda in (0,+inf].
 */
double planckpower(double T);

/**
 * Fits the value between to min and max to range 0 to 1
 */
double fit(double val, double min, double max);

/**
 * Calculate a very close approxumation to the chromaticity (hue) of a
 * blackbody emitter at the given temperature 'T' (in Kelvin).
 * This version is an approximation using 3 cubic splines. Not the
 * exact result, but many times faster to compute, so it's the default.
 * The result is in linear color space(!!) and is modified to have equal 
 * luminance across the entire output range.
 */
void bbspline(double T, double xyz[]);

/**
 * XYZ -> Chromaticity[3] (xyz)
 */
void xyztochroma(double xyz[]);

/**
 * Video luminance (excluding gamma)
 * based on Houdini default:
 *     Luma:Paul Haeberli {0.3086, 0.6094, 0.0820}
 */
double luma(double values[]);

/**
 * Perform tonemapping
 */
void tonemap(double rgb[], double avgloglum, double key, double burn);

/**
 * Returns the log to a given base
 */
double log2(double x, double base);

/**
 * Performs black body color mapping
 * 
 * @param T         - Temperature (K)
 * @param adapt     - tonemapping adaption
 * @param burn      - tonemapping burn
 * @param out_chr   - chromaticity
 * @return          - power/intensity
 * */
double blackbody(double T, double adapt, double burn, double out_chr[]);

#endif //BLACKBODY_H