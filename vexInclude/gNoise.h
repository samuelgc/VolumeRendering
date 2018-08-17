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
 * NAME:	gUtil.h (Gallery Library, VEX)
 *
 * COMMENTS:	This include file contains many useful functions for writing
 *		shaders.
 */

#ifndef __gNoise_h__
#define __gNoise_h__

#include "gUtil.h"

#define gDEFAULT_OCTAVES	12

//
// Generic define for fBm noise
//
#define gNOISE_FUNC(type, nfunc, scaleP)	\
	    float	scale, plimit, blend; \
	    type	nval; \
	    int		octaves; \
	    plimit = 2*blur; \
	    nval = 0; \
	    octaves = 0; \
	    scale = 1; \
	    while (scale > plimit && octaves < maxoct) { \
		nval += scale*(type(nfunc) - .5); \
		scale *= rough; \
		scaleP; \
		octaves++; \
	    } \
	    if (scale > blur) { \
		blend = scale * clamp((scale/blur) - 1, 0, 0); \
		nval += blend * (type(nfunc) - .5); \
	    } \
	    return nval;
	    
//
// gfBmF() computes fractal Brownian motion noise
//	p	= Position to compute noise at
//	amp	= roughness of noise (between 0 and 1)
//	maxoct	= maximum number of octaves of noise to compute (try 10)
//
float
gfBmF(vector p; float rough; int maxoct; float blur)
{
    vector	pp = p;
    gNOISE_FUNC(float, noise(pp), pp *= 2 )
}

float
gfBmFD(vector p; float rough; int maxoct)
{
    return gfBmF(p, rough, maxoct, gFilterWidthV(p));
}

vector
gfBmV(vector p; float rough; int maxoct; float blur)
{
    vector	pp = p;
    gNOISE_FUNC(vector, noise(pp), pp *= 2 )
}

vector
gfBmVD(vector p; float rough; int maxoct)
{
    return gfBmV(p, rough, maxoct, gFilterWidthV(p));
}


//
// 2D versions of fBm noise
//

float
gfBm2F(float u, v; float rough; int maxoct; float blur)
{
    float	x = u, y = v;
    gNOISE_FUNC(float, noise(x, y), x *= 2; y *= 2 )
}

float
gfBm2FD(float u, v; float rough; int maxoct)
{
    return gfBm2F(u, v, rough, maxoct, max(gFilterWidthF(u), gFilterWidthF(v)));
}

vector
gfBm2V(vector u, v; float rough; int maxoct; float blur)
{
    float	x = (float)u;
    float	y = (float)v;
    gNOISE_FUNC(vector, noise(x, y), x *= 2; y *= 2 )
}

vector
gfBm2VD(float u, v; float rough; int maxoct)
{
    return gfBm2V(u, v, rough, maxoct, max(gFilterWidthF(u), gFilterWidthF(v)));
}

//
// Periodic versions of noise functions.  The period must be integral
//
float
gPfBmF(vector p; float rough; int maxoct; int perx, pery, perz; float blur)
{
    vector	pp = p;
    int		px = perx, py = pery, pz = perz;
    gNOISE_FUNC(float, pnoise(pp, px, py, pz), pp *= 2; px *= 2; py *= 2; pz *= 2 )
}

float
gPfBmFD(vector p; float rough; int maxoct; int perx, pery, perz)
{
    return gPfBmF(p, rough, maxoct, perx, pery, perz, gFilterWidthV(p));
}

vector
gPfBmV(vector p; float rough; int maxoct; int perx, pery, perz; float blur)
{
    vector	pp = p;
    int		px = perx, py = pery, pz = perz;
    gNOISE_FUNC(vector, pnoise(pp, px, py, pz), pp *= 2; px*=2; py*=2; pz*=2 )
}

vector
gPfBmVD(vector p; float rough; int maxoct; int perx, pery, perz)
{
    return gPfBmV(p, rough, maxoct, perx, pery, perz, gFilterWidthV(p));
}


//
// 2D versions of periodic fBm noise
//

float
gPfBm2F(float u, v; float rough; int maxoct; int perx, pery; float blur)
{
    float	x = u, y = v;
    int		px = perx, py = pery;
    gNOISE_FUNC(float, pnoise(x, y, px, py), x *= 2; y *= 2; px *= 2; py *= 2 )
}

float
gPfBm2FD(float u, v; float rough; int maxoct; int perx, pery)
{
    return gPfBm2F(u, v, rough, maxoct, perx, pery,
		    max(gFilterWidthF(u), gFilterWidthF(v)));
}

vector
gPfBm2V(vector u, v; float rough; int maxoct; int perx, pery; float blur)
{
    float	x = (float)u;
    float	y = (float)v;
    int		px = perx, py = pery;
    gNOISE_FUNC(vector, pnoise(x, y, px, py), x *= 2; y *= 2; px *= 2; py *= 2 )
}

vector
gPfBm2VD(float u, v; float rough; int maxoct; int perx, pery)
{
    return gPfBm2V(u, v, rough, maxoct, perx, pery,
		max(gFilterWidthF(u), gFilterWidthF(v)));
}

vector
gNoiseOffset3D(vector pp; int nperiodic; vector nfreq, noff; float nrough)
{
    vector	nval;
    if (nperiodic)
    {
	float	px, py, pz;
	vector	infreq;

	infreq = abs(rint(nfreq));
	assign(px, py, pz, infreq);
	nval = gPfBmVD(pp*nfreq+noff, nrough, gDEFAULT_OCTAVES,
				(int)px, (int)py, (int)pz);
    }
    else nval = gfBmVD(pp*nfreq+noff, nrough, gDEFAULT_OCTAVES);

    return nval;
}

#endif
