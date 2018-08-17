/*
 * PROPRIETARY INFORMATION.  This software is proprietary to
 * Side Effects Software Inc., and is not to be reproduced,
 * transmitted, or disclosed in any way without written permission.
 *
 * Produced by:
 *	Mark Elendt
 *	Side Effects Software Inc
 *	477 Richmond Street West
 *	Toronto, Ontario
 *	Canada   M5V 3E7
 *	416-504-9876
 *
 * NAME:	gBump.h ( Gallery Library, VEX)
 *
 * COMMENTS:	Contains various bump routines.  The bump return an "amount" to
 *		displace the surface.
 */

#ifndef __gBump__
#define __gBump__

#include "gNoise.h"
#include "gStamp.h"

float
gBump3DTile(vector p; float groove_width; float blur)
{
    float	scale;
    scale = gPulseTrain(groove_width, p.x, blur);
    scale = max(scale, gPulseTrain(groove_width, p.y, blur));
    scale = max(scale, gPulseTrain(groove_width, p.z, blur));
    return scale;
}

float
gBump3DTileD(vector p; float groove_width)
{
    return gBump3DTile(p, groove_width, gFilterWidthV(p));
}

float
gBump2DTile(float u, v; float groove_width; float blur)
{
    float	scale;
    scale = gPulseTrain(groove_width, u, blur);
    scale = max(scale, gPulseTrain(groove_width, v, blur));
    return scale;
}

float
gBump2DTileD(float u, v; float groove_width)
{
    float	scale;
    scale = gPulseTrain(groove_width, u, gFilterWidthF(u));
    scale = max(scale, gPulseTrain(groove_width, v, gFilterWidthF(v)));
    return scale;
}

float
gBump3DDent(vector p; float rough; int maxoctaves; float blur)
{
    return gfBmF(p, rough, maxoctaves, blur);
}

float
gBump3DDentD(vector p; float rough; int maxoctaves)
{
    return gfBmFD(p, rough, maxoctaves);
}

float
gBump2DDent(float u, v; float rough; int maxoctaves; float blur)
{
    return gfBm2F(u, v, rough, maxoctaves, blur);
}

float
gBump2DDentD(float u, v; float rough; int maxoctaves)
{
    return gfBm2FD(u, v, rough, maxoctaves);
}

float
gBump2DDimple(float u, v; float thresh, rough; int maxoctaves; float blur)
{
    float	nval = gfBm2F(u, v, rough, maxoctaves, blur);
    return clamp(abs(nval), thresh, 1) - thresh;
}

float
gBump2DDimpleD(float u, v; float thresh, rough; int maxoctaves)
{
    float	nval = gfBm2FD(u, v, rough, maxoctaves);
    return clamp(abs(nval), thresh, 1) - thresh;
}

float
gBump3DDimple(vector p; float thresh, rough; int maxoctaves; float blur)
{
    float	nval = gfBmF(p, rough, maxoctaves, blur);
    return clamp(abs(nval), thresh, 1) - thresh;
}

float
gBump3DDimpleD(vector p; float thresh, rough; int maxoctaves)
{
    float	nval = gfBmFD(p, rough, maxoctaves);
    return clamp(abs(nval), thresh, 1) - thresh;
}

float
gBumpMap(float ss, tt; string mapname; float blur)
{
    vector	clr, uclr, vclr;
    clr = texture(mapname, ss-blur, tt-blur,
			   ss+blur, tt-blur,
			   ss+blur, tt+blur,
			   ss-blur, tt+blur, "filter", "gauss");
    uclr = Du(clr);
    vclr = Dv(clr);
    return luminance(uclr) + luminance(vclr);
}

float
gBumpMapD(float ss, tt; string mapname)
{
    return gBumpMap(ss, tt, mapname, max(gFilterWidthF(ss), gFilterWidthF(tt)));
}

//
// Returns the amount to displace a point
//
float
gComputeBump3D(string bumptype;	// Bump style
	  vector pp;		// Shading space position and original position
	  vector bumpfreq;	// Bump frequency
	  vector bumpoff;	// Bump offset
	  float	 groovewidth;	// For tile based textures
	  float	 noiserough;	// For dent based textures
	  string mapname;	// For texture mapping
    )
{
    vector	pshade;
    float	bscale;

    bscale = 0;
    if (bumptype != "none")
    {
	pshade	= pp * bumpfreq + bumpoff;
	if (bumptype == "tile")
	    bscale = gBump3DTileD(pshade, groovewidth);
	else if (bumptype == "dent")
	    bscale = gBump3DDentD(pshade, noiserough, 6);
	else if (bumptype == "dimple")
	    bscale = gBump3DDimpleD(pshade, groovewidth, noiserough, 6);
    }
    return bscale;
}

float
gComputeBump2D(string bumptype;		// Bump style
	      vector	st;		// 2D texture coordinates
	      vector	bumpfreq;	// Bump frequency
	      vector	bumpoff;	// Bump offset
	      float	groovewidth;	// For tile based textures
	      float	noiserough;	// For dent based textures
	      string	mapname;	// For texture mapping
    )
{
    vector	bst;
    float	bscale;

    bscale = 0;
    if (bumptype != "none")
    {
	bst	= st*bumpfreq + bumpoff;
	if (bumptype == "uvtile")
	    bscale = gBump2DTileD(bst.x, bst.y, groovewidth);
	else if (bumptype == "uvdent")
	    bscale = gBump2DDentD(bst.x, bst.y, noiserough, 6);
	else if (bumptype == "uvdimple")
	    bscale = gBump2DDimpleD(bst.x, bst.y, groovewidth, noiserough, 6);
	else if (bumptype == "map" && mapname != "")
	    bscale = gBumpMapD(bst.x, bst.y, mapname);
    }
    return bscale;
}

#define gComputeNV3D(nn, vv, P, I, bumptype,pp, bumpamp, bumpfreq, bumpoff, \
		     groovewidth, noiserough, mapname) \
	vv = -normalize(I); \
	if (bumptype != "none") { \
	    vector pptmp; \
	    float bscale = gComputeBump3D(bumptype, pp, \
		    bumpfreq, bumpoff, groovewidth, noiserough, mapname); \
	    pptmp = P + bumpamp*bscale*normalize(N); \
	    nn = normalize(frontface(computenormal(pptmp), I)); \
	} else nn = normalize(frontface(N, I));

#define gComputeNV2D(nn, vv, P, I, bumptype, st, bumpamp, bumpfreq, bumpoff, \
		     groovewidth, noiserough, mapname) \
	vv = -normalize(I); \
	if (bumptype != "none") { \
	    vector pptmp2D; \
	    float bscale2D = gComputeBump2D(bumptype, st, \
		    bumpfreq, bumpoff, groovewidth, noiserough, mapname); \
	    pptmp2D = P + bumpamp*bscale2D*normalize(N); \
	    nn = normalize(frontface(computenormal(pptmp2D), I)); \
	} else nn = normalize(frontface(N, I));

#endif
