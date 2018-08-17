/*
 * PROPRIETARY INFORMATION.  This software is proprietary to
 * Side Effects Software Inc., and is not to be reproduced,
 * transmitted, or disclosed in any way without written permission.
 *
 * Produced by:
 *	Side Effects Software Inc
 *	477 Richmond Street West
 *	Toronto, Ontario
 *	Canada   M5V 3E7
 *	416-504-9876
 *
 * NAME:	shading.h (VEX)
 *
 * COMMENTS:	This include file contains useful shading functions
 */

#ifndef __shading_h__
#define __shading_h__

#include <math.h>

//
// Defines for the mask value passed to the illumanance() construct.
//	illuminance(vector position, vector axis, float angle, int mask)
//	
#define LIGHT_AMBIENT	bouncemask("ambient")
#define LIGHT_DIFFUSE	bouncemask("diffuse")
#define LIGHT_SPECULAR	bouncemask("reflect")
#define LIGHT_DIFFSPEC	(LIGHT_DIFFUSE|LIGHT_SPECULAR)

//
// This method will scale the angle between dir and axis, making it
// possible to squash or expand environment maps about an axis.
//
vector
shading_computeEnvAngleScale(vector dir; vector axis; float anglescale)
{
    vector	ndir;

    if (anglescale != 1)
    {
	float	zangle = acos(dir.z);
	vector	raxis;
	matrix3	rmat;

	raxis = normalize(cross(axis, dir));
	if (anglescale > 0)
	{
	    zangle /= anglescale;
	    zangle = clamp(zangle, 0.0, PI);
	}
	else
	    zangle = PI;

	rmat = ident();
	rotate(rmat, zangle, raxis);

	ndir = axis * rmat;
    }
    else
	ndir = dir;

    return ndir;
}

// Conductor fresnel functions
vector fresnelcond_nmin(vector r)
{
    return (1-r)/(1+r);
}

vector fresnelcond_nmax(vector r)
{
    return (1+sqrt(r))/(1-sqrt(r));
}

vector fresnelcond_eta(vector r, g)
{
    return fresnelcond_nmin(r)*g + (1-g)*fresnelcond_nmax(r);
}

vector fresnelcond_kappa2(vector r, n)
{
    vector nr = (n+1)*(n+1)*r-(n-1)*(n-1);
    return nr/(1-r);
}

vector fresnelcond_reflectivity(vector n, k)
{
    return ((n-1)*(n-1)+k*k)/((n+1)*(n+1)+k*k);
}

vector fresnelcond_edgetint(vector n, r)
{
    return (fresnelcond_nmax(r)-n)/(fresnelcond_nmax(r)-fresnelcond_nmin(r));
}

// Artistic variant with Reflectivity, Edgetint parameterization
vector fresnelcond_artistic(vector nI, nN, r, g)
{
    float dotNI = clamp(dot(nN, -nI), 0, 0.999);

    vector n = fresnelcond_eta(r, g);
    vector k2 = fresnelcond_kappa2(r, n);

    vector rs_num = n*n + k2 - 2*n*dotNI + dotNI*dotNI;
    vector rs_den = n*n + k2 + 2*n*dotNI + dotNI*dotNI;
    vector rs = rs_num/rs_den;

    vector rp_num = (n*n + k2)*dotNI*dotNI - 2*n*dotNI + 1;
    vector rp_den = (n*n + k2)*dotNI*dotNI + 2*n*dotNI + 1;
    vector rp = rp_num/rp_den;

    return 0.5*(rs+rp);
}

// Physical variant with pysical parameters
// Eta (Refractive Index) and
// Kappa (Extinction Coefficient)
vector fresnelcond_physical(vector nI, nN, eta, kappa)
{
    float dotNI = clamp(dot(nN, -nI), 0, 0.999);

    vector tmp = (eta*eta + kappa*kappa) * dotNI*dotNI;
    vector Rparl2 = (tmp - (2.0f * eta * dotNI) + 1) /
                    (tmp + (2.0f * eta * dotNI) + 1);
    vector tmp_f = eta*eta + kappa*kappa;
    vector Rperp2 = (tmp_f - (2.0f * eta * dotNI) + dotNI*dotNI) /
                    (tmp_f + (2.0f * eta * dotNI) + dotNI*dotNI);
    return 0.5 * (Rparl2 + Rperp2);
}

void thinfresnel(const vector nI, nN; const float eta; export float kr, kt)
{
    float kr1, kt1;
    vector R1, T1;
    fresnel(nI, nN, eta, kr1, kt1, R1, T1);

    float kr2, kt2;
    vector R2, T2;
    fresnel(T1, nN, eta, kr2, kt2);

    float tmp = 1.0/(1 - kr2*kr2);
    kr = kr1 + (kt1*kr2*kt2) * tmp;
    kt = kt1*kt2*tmp;
}

void thinfresnel(const vector nI, nN; const float eta; export float kr, kt; export vector R, T)
{
    thinfresnel(nI, nN, eta, kr, kt);

    R = reflect(nI, nN);
    T = nI;
}

#if defined(VOP_SHADING)
float
dirtmask_dome(vector P;
	 vector nN;
	 int inSID;
	 int	samples;
	 float  maxdist;
	 vector biasdir;
	 float  bias;
	 string scope)
{
    int hit, found;
    float dist;
    vector localdir, dir, hitP, hitN;

    vector u = normalize(set(nN.z, nN.z, -nN.x-nN.y));
    vector v = normalize(cross(nN, u));
    matrix3 space = set(u, v, nN);

    float _maxdist;
    vector _biasdir = biasdir;
    _biasdir = ntransform("space:world", "space:current", _biasdir);
    _biasdir = normalize(select(bias >= 0, -_biasdir, _biasdir));
    float _bias = abs(bias);

    int hits = 0;
    vector2 sample;

    float raybias;
    found = renderstate("renderer:raybias", raybias);
    raybias = select(found, raybias, 1e-3);

    string rengine;
    renderstate("renderer:renderengine", rengine);
    int israytrace = rengine == "raytrace" || rengine == "pbrraytrace";
    int sid = israytrace ? inSID : newsampler();

    for(int i=0; i<samples; i++)
    {
	if (israytrace)
	    nextsample(sid, sample.x, sample.y, "mode", "nextpixel");
	else
	    nextsample(sid, sample.x, sample.y, "mode", "qstrat");

	float phi = 2.0*M_PI*sample.x;
	float theta = asin(sqrt(sample.y));

	localdir.x = cos(theta) * sin(phi);
	localdir.y = cos(theta) * cos(phi);
	localdir.z = sin(theta);
	dir = localdir * space;

	float dotbias = dot(dir, _biasdir);
	float dotNB = dot(nN, _biasdir);
	_biasdir = _biasdir - clamp(dotNB, -1, 0) * nN;
	dir += _bias * clamp(dotbias, 0, 1) * _biasdir;
	dir *= fit(dotbias, -2, 1, 1.0/(1.0+_bias), 1);

	vector testNg = select(dot(nN, Ng) > 0.0, Ng, -Ng);
	if (dot(dir, testNg) < 1e-5)
		continue;

	dist = rayhittest(P, maxdist * dir, hitP, hitN, raybias, "scope", scope);
	hit = dist >= 0.0;

	hits += hit;
    }

    return (float)hits/samples;
}

float
dirtmask_fan(vector P;
	 vector nN;
	 int inSID;
	 float  sensitivity;
	 int	samples;
	 float  maxdist;
	 vector biasdir;
	 float  bias;
	 string scope)
 {
    int hit, found;
    float dist;
    vector localdir, dir, hitP, hitN;

    vector u = normalize(set(nN.z, nN.z, -nN.x-nN.y));
    vector v = normalize(cross(nN, u));
    matrix3 space = set(u, v, nN);

    float max_dirt = 0.0;

    vector2 sample;
    float _maxdist;
    vector _biasdir = normalize(select(bias >= 0, -biasdir, biasdir));
    _biasdir = ntransform("space:world", "space:current", _biasdir);
    float _bias = abs(bias);

    float raybias;
    found = renderstate("renderer:raybias", raybias);
    raybias = select(found, raybias, 1e-3);

    string rengine;
    renderstate("renderer:renderengine", rengine);
    int israytrace = rengine == "raytrace" || rengine == "pbrraytrace";
    int sid = israytrace ? inSID : newsampler();

    for(int i=0; i<samples; i++)
    {
	if (israytrace)
	    nextsample(sid, sample.x, sample.y, "mode", "nextpixel");
	else
	    nextsample(sid, sample.x, sample.y, "mode", "qstrat");

	float phi = 2.0 * M_PI * sample.x;
	float theta = 0.5*M_PI*sensitivity;

	localdir.x = cos(phi);
	localdir.y = sin(phi);
	localdir.z = 0;

	dir = localdir * space;

	float dotNB = dot(nN, _biasdir);
	_biasdir = _biasdir - dotNB * nN;
	float dotbias = dot(dir, normalize(_biasdir));
	dir += _bias * clamp(dotbias, 0, 1) * _biasdir;
	dir *= fit(dotbias, -2, 1, 1.0/(1.0+_bias), 1);

	float len = length(dir);
	dir *= sin(theta);
	dir += fit01(sin(theta), 1, len) * cos(theta) * nN;

	dist = rayhittest(P, maxdist * dir, hitP, hitN, raybias, "scope", scope);

	hit = dist >= 0.0;

	float dirt = clamp(1.0 - dist / (maxdist * length(dir)), 0, 1);
	max_dirt = select(hit, max(dirt, max_dirt), max_dirt);
    }

    return max_dirt;
}

#endif

#endif
