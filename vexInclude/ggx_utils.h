#ifndef __ggx_utils_h__
#define __ggx_utils_h__

#include <math.h>

int SameHemisphere(vector w; vector wp) {
    return (w.z * wp.z) > 0.0f;
}

float sqr(float v) {
    return v*v;
}

float schlick_f(float u) {
    float m = clamp(1.0-u, 0.0, 1.0);
    float m2 = m*m;
    return m2*m2*m;
}

float smith_g(vector w; float alphax, alphay) {
    if (w.z <= 0.0) {
        return 0.0;
    }

    float tanTheta2 = (1 - w.z * w.z) / (w.z * w.z);
    float cosPhi = w.x;
    float sinPhi = w.y;

    float alpha2 = (cosPhi*cosPhi)*(alphax*alphax) + (sinPhi*sinPhi)*(alphay*alphay);
    alpha2 /= cosPhi*cosPhi + sinPhi*sinPhi;

    return 2.0 / (1.0 + sqrt(1.0 + alpha2 * tanTheta2));
}

vector ggx_albedo_estimate(const float cosTheta, rough, ior; const vector F0)
{
    float cosNI = 1.0 - abs(cosTheta);
    float tmp2 = pow(cosNI, lerp(5, 1.9, rough));

    float divide3 = sqr(ior-1) / sqr(ior+1);

    vector facing = F0 * fit01(pow(rough, 2.7), divide3, 0.31*divide3);
    float pow2 = pow(rough, 1-rough);
    vector grazing = (vector) lerp(1, rough * divide3 * 0.9, pow2);
    return select(ior > 1.000001, lerp(facing, grazing, tmp2), 0.0);
}

// http://graphicrants.blogspot.de/2013/08/specular-brdf-reference.html eq. 5
float computeD(vector wm; float alphax, alphay)
{
    return 1.0 / sqr(sqr(wm.x/alphax) + sqr(wm.y/alphay) + sqr(wm.z));
}

// Importance Sampling Microfacet-Based BSDFs using
// the Distribution of Visible Normals - Eric Heitz, Eugene D’Eon
//
// This converts the distribution 'D' to the distribution 'Dwi' of normals
// visible from wi.
float computeDwi(vector wi, wm; float Gi, D)
{
    // Eq. 2
    return Gi * abs(dot(wi, wm)) * D / abs(wi.z);
}

vector2 sample_ggx_slope(const float theta_i, sx, sy)
{
    vector2 slope;

    // special case (normal incidence)
    if(theta_i < 0.0001)
    {
	float r = sqrt(sx/(1-sx));
	float phi = M_TWO_PI * sy;
	slope.x = r * cos(phi);
	slope.y = r * sin(phi);
	return slope;
    }

    // precomputations
    float tan_theta_i = tan(theta_i);
    float a = 1 / (tan_theta_i);
    float G1 = 2 / (1 + sqrt(1.0+1.0/(a*a)));

    // sample slope.x
    float A = 2.0*sx/G1 - 1.0;
    float tmp = 1.0 / (A*A-1.0);
    float B = tan_theta_i;
    float D = sqrt(B*B*tmp*tmp - (A*A-B*B)*tmp);
    float slope_x_1 = B*tmp - D;
    float slope_x_2 = B*tmp + D;
    slope.x = (A < 0 || slope_x_2 > 1.0/tan_theta_i) ? slope_x_1 : slope_x_2;

    // sample slope.y
    float S = select(sy > 0.5, 1.0, -1.0);
    float sy_ = select(sy > 0.5, 2.0*(sy-0.5), 2.0*(0.5-sy));

    float z = (sy_*(sy_*(sy_*0.27385-0.73369)+0.46341)) / (sy_*(sy_*(sy_*0.093073+0.309420)-1.000000)+0.597999);
    slope.y = S * z * sqrt(1.0+slope.x*slope.x);

    return slope;
}

vector sample_ggx(const vector wi; const float alpha_x, alpha_y, sx, sy)
{
    // stretch wi
    vector wi_ = normalize(set(alpha_x * wi.x, alpha_y * wi.y, wi.z));

    // get polar coordinates of wi_
    int normalinc = wi_.z >= 0.99999;
    float theta_ = select(normalinc, 0.0, acos(wi_.z));
    float phi_ = select(normalinc, 0.0, atan2(wi_.y, wi_.x));

    // sample
    vector2 slope = sample_ggx_slope(theta_, sx, sy);

    // rotate
    float tmp = cos(phi_)*slope.x - sin(phi_)*slope.y;
    slope.y = sin(phi_)*slope.x + cos(phi_)*slope.y;
    slope.x = tmp;

    // unstretch
    slope.x = alpha_x * slope.x;
    slope.y = alpha_y * slope.y;

    // compute normal
    return normalize(set(-slope.x, -slope.y, 1.0));
}

#endif
