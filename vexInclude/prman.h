#ifndef __prman_h__
#define __prman_h__

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
 * NAME:	prman.h (VEX)
 *
 * COMMENTS:	This file contains defines and inline functions which are
 *		helpful when porting RenderMan(tm) shaders to VEX code.
 *
 * NOTES:
 *	- The shader message passing functions are implemented in VEX via:
 *		lightsource()	:= int limport(string name, {type}data)
 *		surface()	:= int simport(string name, {type}data)
 *		displacement()	:= int dimport(string name, {type}data)
 *
 *	  Where {type} is one of:
 *		int, float, vector, vector4, matrix3, matrix
 *
 *	  Since atmosphere shading is done last in the pipeline, there are no
 *	  shaders left to read export variables, so there is no VEX
 *	  correspondence to the prman atmosphere() function.
 */

#include <math.h>

// Shader name mappings
#define volume	fog

// Type mappings
#define color	vector
#define	normal	vector
#define point	vector
#define hpoint	vector4

// Storage mappings
#define uniform
#define varying
#define output	export

// Global variable mappings
#define Ci			Cf
#define Oi			Of

// WARNING, users may want to undefine these and pass in parameters for VEX
#define Os			{1,1,1}
#define Cs			{1,1,1}

#undef E	// Undefine the constant E and re-define the variable
#define E			Eye

// Function mappings
#define	faceforward		frontface
#define mix(aa,bb,cc)		lerp(aa,bb,cc)
#define step(min,value)		((value) > (min))
#define smoothstep(aa,bb,cc)	smooth(aa,bb,cc)
#define round(aa)		rint(aa)
#define spline			cspline
#define illuminate(xx,yy,zz)	
#define solar(xx,yy)		{ L = -Lz; }
#define xcomp(aa)		(aa.x)
#define ycomp(aa)		(aa.y)
#define zcomp(aa)		(aa.z)
#define	cellnoise(aa)		random(aa)

//
// WARNING:  To be correct, the magnitude of dPds needs to be scaled by Du(s).
// However, since this is a macro, the computation will be performed on every
// evaluation.  If the variable is accessed more than one time, it would be
// better to assign a local.
//
#define dPdu			(dPds/Du(s))
#define dPdv			(dPdt/Dv(t))
#define calculatenormal		computenormal

float
prman_mod(float aa; float bb)
{
    float	rem = aa % bb;
    if (rem < 0) rem += bb;
    return rem;
}

#define mod(aa, bb)	prman_mod((aa), (bb))

#endif
