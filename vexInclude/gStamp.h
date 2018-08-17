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
 * NAME:	gStamp.h (Gallery Library, VEX)
 *
 * COMMENTS:	This include file contains many stamping functions.
 */

#ifndef __gStamp_h__
#define __gStamp_h__

#include "gUtil.h"

//
// Pulse:  returns the "filtered" version of (x >= edge0 && x < edge1)
//	edge0	= leading edge of the pulse
//	edge1	= trailing edge of the pulse
//	x	= sample point (returns 1 if it's in the pulse
//	blur	= softness of pulse (default 1)
//	filter	= filter kernel name
//
float
gPulse(float edge0, edge1, x; float blur; string filter)
{
    float	x0, x1, e0, e1;

    x0 = x - blur;
    x1 = x + blur;
    e0 = filterstep(edge0, x0, x1, "filter", filter);
    e1 = filterstep(edge1, x0, x1, "filter", filter);
    return e0 - e1;
}

float
gPulseD(float edge0, edge1, x; string filter)
{
    return gPulse(edge0, edge1, x, 0.5*gFilterWidthF(x), filter);
}

float
gPulseTrain(float edge0, x; float blur)
{
    float	x0, x1, e1;
    x0 = x - blur;
    x1 = x + blur;
    e1 = 1 - edge0;
    x0 = edge0*floor(x0) + max(0, frac(x0) - e1);
    x1 = edge0*floor(x1) + max(0, frac(x1) - e1);
    return (x1 - x0) / blur;
}

float
gPulseTrainD(float edge0, x)
{
    return gPulseTrain(edge0, x, 0.5*gFilterWidthF(x));
}

//
// The checkerboard function requires $HH/map/gChecker.rat to exist.  This
// should be a black and white checkerboard picture resolution 64x64.
float
gCheckerBoard(float x, y; float blur0, blur1)
{
    return vector(texture("gBWCheck.rat", x - blur0, y - blur1,
					 x + blur0, y - blur1,
					 x + blur0, y + blur1,
					 x - blur0, y + blur1,
					 "filter", "gauss")).r;
}

float
gCheckerBoardD(float x, y)
{
    return gCheckerBoard(x, y, .5*gFilterWidthF(x), .5*gFilterWidthF(y));
}


//
// Circle Stamp:  Returns the "filtered" version of (px*px + py*py) < 1
//	px	= Sample point (relative to the unit circle)
//	py	= Y coordinate
//	blur	= softness of circle (default 1)
//	filter	= filter kernel name
float
gCircleStamp(float px, py; float blur; string filter)
{
    float	d, dd;
    d = px*px + py*py;
    return 1 - filterstep(1, d-blur, d+blur, "filter", filter);
}

float
gCircleStampD(float px, py; string filter)
{
    float	fx, fy;
    fx = gFilterWidthF(px);
    fy = gFilterWidthF(py);
    return gCircleStamp(px, py, max(fx, fy), filter);
}

//
// Circle Stamp:  Returns the "filtered" version of (px*px + py*py) < 1
//	px	= Sample point (relative to the unit circle)
//	py	= Y coordinate
//	iradius	= inside radius
//	oradius	= outside radius
//	blur	= softness of circle (default 1)
//	filter	= filter kernel name
float
gRingStamp(float px, py; float iradius, oradius; float blur; string filter)
{
    float	d, dd;
    d = px*px + py*py;
    return 1 - gPulse(iradius, oradius, d, blur, filter);
}

float
gRingStampD(float px, py; float iradius, oradius; string filter)
{
    float	fx, fy;
    fx = gFilterWidthF(px);
    fy = gFilterWidthF(py);
    return gRingStamp(px, py, iradius, oradius, max(fx, fy), filter);
}

//
// Sphere Stamp:  Returns the "filtered" version of (length2(p) < 1)
//	p	= Sample point (relative to unit sphere)
//	blur	= filter softness (default 1)
//	filter	= filter kernel name
//
float
gSphereStamp(vector p; float blur; string filter)
{
    float	d;

    d = length2(p);
    return 1 - filterstep(1.0, d-blur, d+blur, "filter", filter);
}

float
gSphereStampD(vector p; string filter)
{
    return gSphereStamp(p, gFilterWidthV(p), filter);
}

//
// Shell Stamp:  Returns the "filtered" version of
//		 (length(p) < oradius) && (length(p) > iradius)
//	p	= Sample point (relative to unit sphere)
//	blur	= filter softness (default 1)
//	filter	= filter kernel name
//
float
gShellStamp(vector p; float iradius, oradius; float blur; string filter)
{
    return gPulse(iradius, oradius, length(p), blur, filter);
}

float
gShellStampD(vector p; float iradius, oradius; string filter)
{
    return gShellStamp(p, iradius, oradius, gFilterWidthV(p), filter);
}

// 
// Box stamp to [-1,-1]-[1,1]:
//	Returns whether the sample point is in the unit box
//	x	= X coordinate
//	y	= Y coordinate
//	blur	= filter softness (default 1)
//	filter	= filter kernel name
//
float
gBox2DStamp(float x, y; float blur; string filter)
{
    float	dx;
    dx  = gPulse(-1, 1, x, blur, filter);
    dx *= gPulse(-1, 1, y, blur, filter);
    return dx;
}

float
gBox2DStampD(float x, y; string filter)
{
    float	dx;
    dx  = gPulseD(-1, 1, x, filter);
    dx *= gPulseD(-1, 1, y, filter);
    return dx;
}

// 
// Box stamp to [-1,-1,-1]-[1,1,1]:
//	Returns whether the sample point is in the unit box
//	x	= X coordinate
//	y	= Y coordinate
//	blur	= filter softness (default 1)
//	filter	= filter kernel name
//
float
gBox3DStamp(vector p; float blur; string filter)
{
    float	dx;
    dx  = gPulse(-1, 1, p.x, blur, filter);
    dx *= gPulse(-1, 1, p.y, blur, filter);
    dx *= gPulse(-1, 1, p.z, blur, filter);
    return dx;
}

float
gBox3DStampD(vector p; string filter)
{
    float	dx;
    dx  = gPulseD(-1, 1, p.x, filter);
    dx *= gPulseD(-1, 1, p.y, filter);
    dx *= gPulseD(-1, 1, p.z, filter);
    return dx;
}

#endif
