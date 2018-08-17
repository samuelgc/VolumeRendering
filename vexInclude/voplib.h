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
 * NAME:	voplib.h ( VOP Library Functions, VEX )
 *
 * COMMENTS:	Support functions for VOPS.  Most of these library functions
 *		are shading context only.
 */

#ifndef __voplib__
#define __voplib__

#if !defined(VOP_SHADING) && !defined(VOP_OP)
    #define VOP_OP
#endif

#define VOP_MIN_FILTER_SIZE	1e-6
#define VOP_PERIOD_FADE		(1.0/8.0)

#if defined(VOP_SHADING)
    #define FILTERSTEP(e0, x0, x1, filter) filterstep(e0,x0,x1,"filter",filter)
#else
    #define FILTERSTEP(e0, x0, x1, filter) clamp(fit(e0, x0, x1, 1, 0), 0, 1)
#endif

#if defined(VOP_SHADING) || defined(VOP_COP2)
    #define AREA(pp)	area(pp)
#else
    #define AREA(pp)	VOP_MIN_FILTER_SIZE
#endif

#define VOP_FADE_PERIODIC(v1, v0, average, fwidth)	\
	    lerp((v1-v0)/fwidth, average, clamp(fwidth*VOP_PERIOD_FADE, 0, 1));

#include <math.h>
#include <shading.h>

#define vop_abs		abs
#define vop_sqrt	sqrt
#define vop_ceil	ceil
#define vop_floor	floor
#define vop_radians	radians
#define vop_degrees	degrees
#define vop_exp		exp
#define vop_frac	frac
#define vop_avg		avg
#define vop_mincomp	min
#define vop_maxcomp	max
#define vop_rint	rint
#define vop_acos	acos
#define vop_asin	asin
#define vop_atan	atan
#define vop_cosh	cosh
#define vop_sinh	sinh
#define vop_tanh	tanh
#define vop_cos		cos
#define vop_sin		sin
#define vop_tan		tan
#define vop_fit		fit
#define vop_log		log
#define vop_smooth	smooth
#define vop_getcomp	getcomp
#define vop_assign	assign
#define vop_vnoise	vnoise
#define vop_frontface	frontface
#define vop_dot		dot
#define vop_cross	cross
#define vop_getrayweight()	getrayweight()
#define vop_luminance(c)	luminance(c)

#define vop_floattovec2(x,y)	set(x,y)
#define vop_floattovec(x,y,z)	set(x,y,z)
#define vop_vec2tofloat(v, x,y)	assign(x,y, v)
#define vop_vectofloat(v, x,y,z)	assign(x,y,z, v)
#define vop_hvectofloat(v, x,y,z,w)	assign(x,y,z,w, v)
#define vop_floattomatrix(a00,a01,a02,a03, a10,a11,a12,a13, \
			  a20,a21,a22,a23, a30,a31,a32,a33) \
		set(a00,a01,a02,a03, a10,a11,a12,a13, \
		    a20,a21,a22,a23, a30,a31,a32,a33)

#define vop_mod(x,y)		((x) % (y))
#define vop_pow(x,y)		pow(x, y)
#define vop_is_even(x)		((int(floor(x)) & 1) == 0)
#define vop_is_odd(x)		((int(floor(x)) & 1) == 1)

/// Random values in a Guassian/Normal distribution using the Box-Muller
/// transform.
float
vop_grandom(float u0, u1; float sigma)
{
    // Given two independent random variables u0 and u1, use the Box-Muller
    // transformation to generate a Gaussian random distribution with the given
    // standard deviation.
    return cos(u0 * (M_PI*2)) * sqrt(-2 * log(u1)) * sigma;
}
float
vop_grandom(int seed; float sigma)
{
    float	u0 = random(seed);
    float	u1 = random(seed*3929 + 48311);
    return vop_grandom(u0, u1, sigma);
}

float
vop_gnrandom(float sigma)
{
    float	u0 = nrandom();
    float	u1 = nrandom();
    return vop_grandom(u0, u1, sigma);
}

vector
vop_grandom(int seed; float sigma)
{
    // Random vector for the given seed, fit to a normal distribution with the
    // given standard deviation.
    float	x = float(vop_grandom(seed, sigma));
    float	y = float(vop_grandom(seed*17 + 97, sigma));
    float	z = float(vop_grandom(seed*29 + 383, sigma));
    return set(x, y, z);
}

vector
vop_gnrandom(float sigma)
{
    // Random vector fit to a normal distribution with the given standard
    // deviation.
    float	x = float(vop_gnrandom(sigma));
    float	y = float(vop_gnrandom(sigma));
    float	z = float(vop_gnrandom(sigma));
    return set(x, y, z);
}

void
vop_bindST(float news, ss; int isSConnected;
	   float newt, tt; int isTConnected)
{
#if defined(VOP_SOP) || defined(VOP_POP)
    news = isSConnected ? ss : P.x;
    newt = isTConnected ? tt : P.y;
#elif defined(VOP_COP) || defined(VOP_COP2)
    news = isSConnected ? ss : X;
    newt = isTConnected ? tt : Y;
#elif defined(VOP_SHADING)
    news = isSConnected ? ss : s;
    newt = isTConnected ? tt : t;
#else
    news = isSConnected ? ss : 0;
    newt = isTConnected ? tt : 0;
#endif
}

void
vop_bindS(float news, ss; int isSConnected)
{
#if defined(VOP_SOP) || defined(VOP_POP)
    news = isSConnected ? ss : P.x;
#elif defined(VOP_COP) || defined(VOP_COP2)
    news = isSConnected ? ss : X;
#elif defined(VOP_SHADING)
    news = isSConnected ? ss : s;
#else
    news = isSConnected ? ss : 0;
#endif
}

void
vop_bindU(export float newu; const float u; int isUVConnected)
{
#if defined(VOP_SOP) || defined(VOP_POP)
    newu = isUVConnected ? u : P.x;
#elif defined(VOP_COP) || defined(VOP_COP2)
    newu = isUVConnected ? u : X;
#elif defined(VOP_SHADING)
    newu = isUVConnected ? u : s;
#else
    newu = isUVConnected ? u : 0;
#endif
}

void
vop_bindUV(export vector2 newuv; const vector2 uv; int isUVConnected)
{
#if defined(VOP_SOP) || defined(VOP_POP)
    newuv = isUVConnected ? uv : set(P.x, P.y);
#elif defined(VOP_COP) || defined(VOP_COP2)
    newuv = isUVConnected ? uv : set(X, Y);
#elif defined(VOP_SHADING)
    newuv = isUVConnected ? uv : set(s, t);
#else
    newuv = isUVConnected ? uv : set(0, 0);
#endif
}

void
vop_bindUV(export vector newuv; const vector uv; int isUVConnected)
{
#if defined(VOP_SOP) || defined(VOP_POP)
    newuv = isUVConnected ? uv : P;
#elif defined(VOP_COP) || defined(VOP_COP2)
    newuv = isUVConnected ? uv : set(X, Y, 0);
#elif defined(VOP_SHADING)
    newuv = isUVConnected ? uv : set(s, t, 0);
#else
    newuv = isUVConnected ? uv : set(0, 0, 0);
#endif
}

///////////////////////////////////////////////////////////////
// gain function from Perlin's section in Texture and Modeling
///////////////////////////////////////////////////////////////
float
vop_bias(float base, bias) 
{
    float val;
    if (base <= 0)
	val = 0;
    else if (base >= 1)
	val = 1;
    else
	val = bias / (((1.0 / base) - 2) * (1 - bias) + 1);
    return val;
}

float
vop_gain(float base, gain) 
{
    float	val;
    if (base < 0.5)
	val = vop_bias(2*base, gain)*.5;
    else
	val = 1-vop_bias(2*(1-base), gain)*.5;
    return val;
}

vector
vop_colorLinearTransform(vector C;
                         string fromSpace;)
{
    vector Ct = C;

    if (fromSpace == "srgb")
    {
        for (int i = 0; i < 3; i++)
        {
            if (C[i] <= 0.04045) Ct[i] = C[i] / 12.92;
            else Ct[i] = pow((C[i] + 0.055) / 1.055, 2.4);
        }
    }

    else if (fromSpace == "linear")
    {
        for (int i = 0; i < 3; i++)
        {
            if (C[i] <= 0.0031308) Ct[i] = C[i] * 12.92;
            else Ct[i] = 1.055 * pow(C[i], 1.0/2.4) - 0.055;
        }
    }

    return Ct;
}

vector
vop_colormix(vector c1, c2; float bias; int adjust)
{
    vector	clr;
    if (adjust == 3)
	clr = cspline(bias, c1, c1, c2, c2);
    else if (adjust == 2)
	clr = lerp(c1, c2, float(smooth(0, 1, bias)));
    else if (adjust == 1)
	clr = lerp(c1, c2, float(clamp(bias, 0, 1)));
    else
	clr = lerp(c1, c2, bias);
    return clr;
}

void
vop_composite(string operation;
	      vector C; float Ca;		// Resulting color/alpha
	      vector A; float Aa;		// Src color (overlay)
	      vector B; float Ba;		// Dest color (background)
	     )
{
    if (operation == "AoverB")
    {
	C  = A  + (1-Aa)*B;
	Ca = Aa + (1-Aa)*Ba;
    }
    else if (operation == "AinsideB")	// A inside B
    {
	C  = A*Ba;
	Ca = Aa*Ba;
    }
    else if (operation == "AoutsideB")
    {
	C  = A*(1-Ba);
	Ca = Aa*(1-Ba);
    }
    else if (operation == "AatopB")
    {
	C  = A*Ba + B*(1-Aa);
	Ca = Ba;
    }
    else if (operation == "AxorB")
    {
	C  = A*(1-Ba) + B*(1-Aa);
	Ca = Aa + Ba - 2*(Aa*Ba);
    }
    else if (operation == "A")
    {
	C  = A;
	Ca = Aa;
    }
    else if (operation == "B")
    {
	C  = B;
	Ca = Ba;
    }
    else if (operation == "clear")
    {
	C = 0;
	Ca = 0;
    }
}

float
vop_ptlined(vector P1, P2, Q)
{
    return ptlined(P1, P2, Q);
}


float
vop_FilterWidth(float x)
{
#if defined(VOP_SHADING) || defined(VOP_COP2)
    float du = Du(x);
    float dv = Dv(x);
    return max(sqrt(du*du + dv*dv), VOP_MIN_FILTER_SIZE);
#else
    // There's really no way to determine what the filter size should be...
    return 0.005;
#endif
}

// Calculate the approxiate filter angle for a given direction vector using
// derivatives.
float
vop_FilterAngle(vector dir)
{
    vector	du, dv;
    vector	a, b, c, d;
    vector	ac, bd;
    float	angle;

    du = 0.5*Du(dir);
    dv = 0.5*Dv(dir);

    a = dir-du-dv;
    b = dir+du-dv;
    c = dir+du+dv;
    d = dir-du+dv;

    ac = cross(a, c) / sqrt(length2(a)*length2(c));
    bd = cross(b, d) / sqrt(length2(b)*length2(d));
    angle = 0.5 * length(cross(ac, bd));

    return acos(1 - (0.5 / PI) * angle);
}

float
vop_Pulse(float edge0, edge1, x, fwidth; string filter)
{
    float	x0, x1;

    x0 = x  - fwidth*.5;
    x1 = x0 + fwidth;
    return max(0, (min(x1, edge1)-max(x0, edge0))/fwidth);
}

float
vop_FilteredSin(float x, fwidth)
{
    float	x0, x1;

    x0 = x - fwidth * .5;
    x1 = x + fwidth;
    return (-M_SQRT1_2)*VOP_FADE_PERIODIC(cos(x1), cos(x0), 0, fwidth);
}

float
vop_FilteredSinD(float x)
{
    return vop_FilteredSin(x, vop_FilterWidth(x));
}

float
vop_FilteredCos(float x, fwidth)
{
    float	x0, x1;

    x0 = x - fwidth * .5;
    x1 = x + fwidth;
    return M_SQRT1_2*VOP_FADE_PERIODIC(sin(x1), sin(x0), 0, fwidth);
}

float
vop_FilteredCosD(float x)
{
    return vop_FilteredCos(x, vop_FilterWidth(x));
}

float
vop_PulseD(float edge0, edge1, x; string filter)
{
    return vop_Pulse(edge0, edge1, x, vop_FilterWidth(x), filter);
}

float
vop_PulseTrain(float edge0, x, fwidth)
{
    float	x0, x1;

    x0 = x  - fwidth*.5;
    x1 = x0 + fwidth;
    x0 = edge0*floor(x0) + min(edge0, frac(x0));
    x1 = edge0*floor(x1) + min(edge0, frac(x1));
    return VOP_FADE_PERIODIC(x1, x0, edge0, fwidth);
}

float
vop_PulseTrainD(float edge0, x)
{
    return vop_PulseTrain(edge0, x, vop_FilterWidth(x));
}

float
vop_RampTrain(float x, fwidth)
{
    float	x0, x1, f;
    x0 = x - fwidth*.5;
    x1 = x0 + fwidth;
    f = frac(x0); x0 = float(floor(x0)) + f*f;
    f = frac(x1); x1 = float(floor(x1)) + f*f;
    return .5*VOP_FADE_PERIODIC(x1, x0, .5, fwidth);
}

float
vop_RampTrainD(float x)
{
    return vop_RampTrain(x, vop_FilterWidth(x));
}

float
vop_IntegrateTent(float x)
{
    float	f;
    f = frac(x);
    if (f > .5) f = f*(2 - f) - .5;
    else	f = f*f;
    return .5*floor(x) + f;
}

float
vop_TentTrain(float x, fwidth)
{
    float	x0, x1;
    x0 = x - fwidth*.5;
    x1 = x0 + fwidth;
    x0 = vop_IntegrateTent(x0);
    x1 = vop_IntegrateTent(x1);
    return VOP_FADE_PERIODIC(x1, x0, .5, fwidth);
}

float
vop_TentTrainD(float x)
{
    return vop_RampTrain(x, vop_FilterWidth(x));
}

float
vop_DotStamp(float px, py, fwidth; string filter)
{
    float	d;
    d = px*px + py*py;
    return FILTERSTEP(1, d-fwidth, d+fwidth, filter);
}

float
vop_BumpStamp(float px, py, fwidth; string filter)
{
    float d  = px*px + py*py;
    return 1-smooth(0.0, 1.0, d);
}

float
vop_RoundCosStamp(float px, py, fwidth; float hexness, sides, power;
		    string filter)
{
    float	ss, tt;
    float	x0, x1;

    ss = atan(py, px);
    tt = px*px + py*py + hexness*pow(abs(1-cos(ss*(sides))), power);
    x0 = tt - fwidth*.5;
    x1 = x0 + fwidth;
    return FILTERSTEP(1, x0, x1, filter);
}

float
vop_RoundSinStamp(float px, py, fwidth; float hexness, sides, power;
		    string filter)
{
    float	ss, tt;
    float	x0, x1;

    ss = atan(py, px) - M_PI;
    tt = px*px + py*py + hexness*pow(abs(1-sin(ss*(sides))), power);
    x0 = tt - fwidth*.5;
    x1 = x0 + fwidth;
    return FILTERSTEP(1, x0, x1, filter);
}


float
vop_RingStamp(float px, py, iradius, oradius, fwidth; string filter)
{
    float	d;
    d = px*px + py*py;
    return vop_Pulse(iradius, oradius, d, fwidth, filter);
}

float
vop_BoxStamp(float px, py, fu, fv; string filter)
{
    float	dx;
    dx  = vop_Pulse(-1, 1, px, fu, filter);
    dx *= vop_Pulse(-1, 1, py, fv, filter);
    return dx;
}

//
// Tile generation
//  Input:	u, v	- The position
//		fx, fy	- Frequency
//		ox, oy	- Offset
//		jitter	- Jitter for each row
//  Return:	id	- Random ID for the tile
//  Output:	u, v,	- The coordinates within the tile
//		du, dv	- Derivative information
//  Unchanged:	fx, fy, ox, oy, jitter
//		
int
vop_TileGen(float u, v; float fx, fy, ox, oy, stagger, jitter)
{
    int		row;

    v = v*fy - oy;
    row = floor(v);
    u = u*fx - ox - stagger*row + jitter*random(row);
    return floor(u) + row*1984;
}

int
vop_HexTileGen(float u, v; float fx, fy, ox, oy)
{
    float	left, right;
    int		row, col;

    v   = v*fy - oy;
    row = floor(v);
    v  = frac(v);

    u = u*fx + ox;
    if (row & 1)
	u += .5;

    col = floor(u);
    u = frac(u);

    if (v > .5)
    {
	right = v - .5;
	left  = 1.5 - v;
	if (u > left)
	{
	    if (!(row & 1)) col++;
	    row++;
	    v -= 1;
	    u = fit(u, left, 1, 0, .5);
	}
	else if (u < right)
	{
	    if (row & 1) col--;
	    row++;
	    v -= 1;
	    u = fit(u, 0, right, 0.5, 1);
	}
	else u = fit(u, right, left, 0, 1);
    }
    v = (v*2+1)/3;
    return row * 938 + col;
}

//
// This function generates Worley noise with a border.  Parameters are:
//	Input:	ss = S coord
//		tt = T coord
//		jx = Jitter in S direction
//		jy = Jitter in T direction
//		bwidth = Border width
//		bsoft = Border softness
//	Output:
//		Return code = 0 in border area, 1 in cell interior
//		centerx = X Center of closest node
//		centery = Y Center of closest node
//		seed = A random seed for the closest node
float
vop_aaCell2D(float ss, tt, jx, jy, bwidth, bsoft, centerx, centery;
		float f1, f2; int seed)
{
    float	p2x, p2y, x0;
    float	blur;
    vector	vp1, vp2;

    vnoise(ss, tt, jx, jy, seed, f1, f2, centerx, centery, p2x, p2y);

    vp1 = set(centerx, centery, 0);
    vp2 = set(p2x, p2y, 0);
    x0 = (f2-f1)*(f1+f2) / max(distance(vp1, vp2), VOP_MIN_FILTER_SIZE);

    blur = max(vop_FilterWidth(ss), vop_FilterWidth(tt)) * (1 + bsoft);

    return FILTERSTEP(bwidth, x0-blur, x0+blur, "gauss");
}


//
// Stamping simple shape patterns.  Some patterns take auxilliary data.  These
// are:
//	ring	- aux data specifies the ring width (0 to 1)
//	hex	- aux specifies the "hexness" (0 to 1)
//
float
vop_StampPattern(int   layers;		// Number of layers
		 float dotsize;		// Base dotsize
		 float dj;		// Dotsize jitter
		 float softness;	// Softness
		 float px, jx;		// X position and jitter
		 float py, jy;		// Y position and jitter
		 float aux, jaux;	// Auxilliary data
		 string spottype;	// Type of spot
		 string filter;		// Filter type
		 float	floatseed;	// Return a random seed
	 )
{
    float	fwidth;
    float	lrandom;
    float	cx, cy, loff, dsize;
    float	fx, fy;
    float	djx, djy;
    float	du, dv;
    float	result;
    float	irad, idot;
    vector	pp;
    int		i;

    du = vop_FilterWidth(px) * softness;
    dv = vop_FilterWidth(py) * softness;
    fwidth = max(du, dv);
    result = 0;
    loff = 0;
    for (i = 0; i < layers; i++, loff += .5)
    {
	lrandom = random(i);
	pp = set(floor(px+loff), floor(py+loff),
		    1000.0*lrandom+500) + {.5, .5, .5};
	floatseed = random(pp);

	pp = vector(random(pp)) - .5;
	dsize = dotsize * (1-dj*floatseed);
	idot = 1/dsize;

	djy  = idot*(1 - dsize);
	djx  = (2*clamp(jx, 0, 1)) * djy;
	djy *= (2*clamp(jy, 0, 1));

	cx = pp.x * djx;
	cy = pp.y * djy;
	fx = cx - 2*(frac(px+loff) - .5)*idot;
	fy = cy - 2*(frac(py+loff) - .5)*idot;

	if (spottype == "ring")
	{
	    irad = 1-clamp(aux + (pp.z - .5)*jaux, 0, 1);
	    result = vop_RingStamp(fx, fy, irad, 1, fwidth, filter);
	}
	else if (spottype == "box")
	{
	    result = vop_BoxStamp(fx, fy, du, dv, filter);
	}
	else if (spottype == "hex")
	{
	    result = 1-vop_RoundCosStamp(fx, fy, fwidth, aux, 6, jaux, filter);
	}
	else if (spottype == "star")
	{
	    result = 1-vop_RoundSinStamp(fx, fy, fwidth, aux, 5, jaux, filter);
	}
	else if (spottype == "bump")
	{
	    result = vop_BumpStamp(fx, fy, fwidth, filter);
	}
	else
	{
	    result = 1-vop_DotStamp(fx, fy, fwidth, filter);
	}
    }
    floatseed = (floatseed - 0.5)*32000;
    return result;
}

float
vop_RipplePattern(float x, y, decay, toff)
{
    float	d;

    d = sqrt(x*x + y*y);
    return sin(d-toff) * exp(d*(-decay));
}

// Implied parameters pos, rough, and maxoctaves
#define VOP_FBMNOISE_FUNC(type)	\
	    float	amp, scale, plimit, blend; \
	    int		octaves; \
	    plimit = 2*blur; \
	    nval = 0; octaves = 0; scale = 1; amp = 1; \
	    while (scale > plimit && octaves < maxoctaves) { \
		if (noisetype == "xnoise") { \
		    nval += amp *(type(xnoise(pp)) - 0.5); \
		} else { \
		    nval += amp *(type(noise(pp)) - 0.5); \
		} \
		amp *= rough; \
		scale *= 0.5; \
		pp *= 2; \
		octaves++; \
	    } \
	    if (scale > blur) { \
		blend = amp * clamp(scale/blur - 1, 0, 1); \
		if (noisetype == "xnoise") { \
		    nval += blend*(type(xnoise(pp)) - 0.5); \
		} else { \
		    nval += blend*(type(noise(pp)) - 0.5); \
		} \
	    }

vector
vop_FlowNoiseGradVV(vector pos; float flow; float delta)
{
    vector	result;

    result.x = flownoise(set(pos.x+delta,pos.y,pos.z), flow)
	     - flownoise(set(pos.x-delta,pos.y,pos.z), flow);

    result.y = flownoise(set(pos.x,pos.y+delta,pos.z), flow)
	     - flownoise(set(pos.x,pos.y-delta,pos.z), flow);

    result.z = flownoise(set(pos.x,pos.y,pos.z+delta), flow)
	     - flownoise(set(pos.x,pos.y,pos.z-delta), flow);

    result /= delta * 2;

    return result;
}

vector4
vop_FlowNoiseGradPP(vector4 pos; float flow; float delta)
{
    vector4	result;

    result.x = flownoise(set(pos.x+delta,pos.y,pos.z,pos.w), flow)
	     - flownoise(set(pos.x-delta,pos.y,pos.z,pos.w), flow);

    result.y = flownoise(set(pos.x,pos.y+delta,pos.z,pos.w), flow)
	     - flownoise(set(pos.x,pos.y-delta,pos.z,pos.w), flow);

    result.z = flownoise(set(pos.x,pos.y,pos.z+delta,pos.w), flow)
	     - flownoise(set(pos.x,pos.y,pos.z-delta,pos.w), flow);

    result.w = flownoise(set(pos.x,pos.y,pos.z,pos.w+delta), flow)
	     - flownoise(set(pos.x,pos.y,pos.z,pos.w-delta), flow);

    result /= delta * 2;

    return result;
}

vector
vop_FlowNoiseGrad3V(vector pos; float flow; float delta)
{
    matrix3	result;
    vector	final;
    vector	dx, dy, dz;

    dx = flownoise(set(pos.x+delta,pos.y,pos.z), flow)
	 - flownoise(set(pos.x-delta,pos.y,pos.z), flow);

    dy = flownoise(set(pos.x,pos.y+delta,pos.z), flow)
	 - flownoise(set(pos.x,pos.y-delta,pos.z), flow);

    dz = flownoise(set(pos.x,pos.y,pos.z+delta), flow)
	 - flownoise(set(pos.x,pos.y,pos.z-delta), flow);

    result = set(dx.x, dx.y, dx.z,
	       dy.x, dy.y, dy.z,
	       dz.x, dz.y, dz.z);

    result /= delta * 2;

    // The correct gradient is a matrix, but we can't step along
    // that so we approximate by stepping along the average of
    // the three possible stepping directions.
    final = 0.577735;
    final *= result;

    return final;
}

vector4
vop_FlowNoiseGrad4P(vector4 pos; float flow; float delta)
{
    matrix	result;
    vector4	final;
    vector	dx, dy, dz, dw;

    dx = flownoise(set(pos.x+delta,pos.y,pos.z,pos.w), flow)
	 - flownoise(set(pos.x-delta,pos.y,pos.z,pos.w), flow);

    dy = flownoise(set(pos.x,pos.y+delta,pos.z,pos.w), flow)
	 - flownoise(set(pos.x,pos.y-delta,pos.z,pos.w), flow);

    dz = flownoise(set(pos.x,pos.y,pos.z+delta,pos.w), flow)
	 - flownoise(set(pos.x,pos.y,pos.z-delta,pos.w), flow);

    dw = flownoise(set(pos.x,pos.y,pos.z,pos.w+delta), flow)
	 - flownoise(set(pos.x,pos.y,pos.z,pos.w-delta), flow);

    result = set(dx.x, dx.y, dx.z, 0,
	       dy.x, dy.y, dy.z, 0,
	       dz.x, dz.y, dz.z, 0,
	       dw.x, dw.y, dw.z, 0);

    result /= delta * 2;

    // The correct gradient is a matrix, but we can't step along
    // that so we approximate by stepping along the average of
    // the three possible stepping directions.
    final = 0.5;
    final *= result;

    return final;
}

float
vop_fbmlength_float(float val)
{
    return val;
}

float
vop_fbmlength_vector(vector val)
{
    return length(val);
}

float
vop_fbmlength_vector4(vector4 val)
{
    return length(val);
}

// Implied parameters pos, rough, maxoctaves, flow, flowrate, advect
#define VOP_FBMFLOWNOISE_FUNC(type, gradtype)	\
	    float	amp, scale, plimit, blend, M, flowv; \
	    int		octaves; \
	    type	namount; \
	    plimit = 2*blur; \
	    M = advect; flowv = flow; \
	    nval = 0; octaves = 0; scale = 1; amp = 1; \
	    while (scale > plimit && octaves < maxoctaves) { \
		namount = amp * (type(flownoise(pp, flowv)) - 0.5); \
		nval += namount; \
		amp *= rough; \
		flowv *= flowrate; \
		scale *= 0.5; \
		if (M != 0.0) \
		{ \
		    pp -= M * vop_fbmlength_##type(namount) * gradtype(pp, flow, 0.01); \
		} \
		pp *= 2; \
		M *= advect; \
		octaves++; \
	    } \
	    if (scale > blur) { \
		blend = amp * clamp(scale/blur - 1, 0, 1); \
		nval += blend*(type(flownoise(pp, flowv)) - 0.5); \
	    }

float
vop_fbmNoiseFF(float pos; float rough; int maxoctaves; string noisetype)
{
    float	blur = vop_FilterWidth(pos);
    float	nval, pp = pos;
    VOP_FBMNOISE_FUNC( float )
    return nval;
}

float
vop_fbmNoiseFV(vector pos; float rough; int maxoctaves; string noisetype)
{
    float	blur = sqrt(AREA(pos));
    vector	pp = pos;
    float	nval;
    VOP_FBMNOISE_FUNC( float )
    return nval;
}

float
vop_fbmNoiseFP(vector4 pos; float rough; int maxoctaves; string noisetype)
{
    float	blur = sqrt(AREA((vector)pos));
    vector4	pp = pos;
    float	nval;
    VOP_FBMNOISE_FUNC( float )
    return nval;
}

vector
vop_fbmNoiseVF(float pos; float rough; int maxoctaves; string noisetype)
{
    float	blur = vop_FilterWidth(pos);
    float	pp = pos;
    vector	nval;
    VOP_FBMNOISE_FUNC( vector )
    return nval;
}

vector
vop_fbmNoiseVV(vector pos; float rough; int maxoctaves; string noisetype)
{
    float	blur = sqrt(AREA(pos));
    vector	pp = pos;
    vector	nval;
    VOP_FBMNOISE_FUNC( vector )
    return nval;
}

vector
vop_fbmNoiseVP(vector4 pos; float rough; int maxoctaves; string noisetype)
{
    float	blur = sqrt(AREA((vector)pos));
    vector4	pp = pos;
    vector	nval;
    VOP_FBMNOISE_FUNC( vector )
    return nval;
}

// Flow noise with turbulence & antialiased.

float
vop_fbmFlowNoiseFV(vector pos; float rough; int maxoctaves; float flow, flowrate, advect)
{
    float	blur = sqrt(AREA(pos));
    vector	pp = pos;
    float	nval;
    VOP_FBMFLOWNOISE_FUNC( float, vop_FlowNoiseGradVV )
    return nval;
}

float
vop_fbmFlowNoiseFP(vector4 pos; float rough; int maxoctaves; float flow, flowrate, advect)
{
    float	blur = sqrt(AREA((vector)pos));
    vector4	pp = pos;
    float	nval;
    VOP_FBMFLOWNOISE_FUNC( float, vop_FlowNoiseGradPP )
    return nval;
}

vector
vop_fbmFlowNoiseVV(vector pos; float rough; int maxoctaves; float flow, flowrate, advect)
{
    float	blur = sqrt(AREA(pos));
    vector	pp = pos;
    vector	nval;
    VOP_FBMFLOWNOISE_FUNC( vector, vop_FlowNoiseGrad3V )
    return nval;
}

vector
vop_fbmFlowNoiseVP(vector4 pos; float rough; int maxoctaves; float flow, flowrate, advect)
{
    float	blur = sqrt(AREA((vector)pos));
    vector4	pp = pos;
    vector	nval;
    VOP_FBMFLOWNOISE_FUNC( vector, vop_FlowNoiseGrad4P )
    return nval;
}

#define VOP_DAMPENFBM_FLOAT \
        while (scale > plimit && octave < maxoctaves) { \
        fbm += scale * fit(noise(pp),0.3,0.7,-1.0,1.0); \
        scale *= lacun; \
        pp /= lacun; \
        octave ++; \
        } \
        scale *= lacun; \
        fbm += scale * fit(noise(pp),0.3,0.7,-1.0,1.0); \
        while (scale > aa) { \
        fbm += scale * fit(noise(pp),0.3,0.7,-1.0,1.0); \
        pp *= lacun; \
        scale *= lacun; \
	    } \

float
vop_dampenFbmFF(float pos, freq, offset; float lacun, amp; int maxoctaves)
{
    float pp = pos * freq + offset;
    float fw = vop_FilterWidth(pp);
    float aa = vop_FilterWidth(pos);
    float plimit = fw/lacun;
    float scale = 1;
    int   octave = 0;
    float fbm = 0;
    VOP_DAMPENFBM_FLOAT
    fbm = amp * fit(fbm, -1.25, 1.25, 0, 1);
    return fbm;
}

float
vop_dampenFbmFV(vector pos, freq, offset; float lacun, amp; int maxoctaves)
{
    vector pp = pos * freq + offset;
    float  fw = sqrt(AREA(pp));
    float  aa = sqrt(AREA(pos));
    float  plimit = fw/lacun;
    float  scale = 1;
    int    octave = 0;
    float  fbm = 0;
    VOP_DAMPENFBM_FLOAT
    fbm = amp * fit(fbm, -1.25, 1.25, 0, 1);
    return fbm;
}

float
vop_dampenFbmFP(vector4 pos, freq, offset; float lacun, amp; int maxoctaves)
{
    vector4 pp = pos * freq + offset;
    float   fw = sqrt(AREA((vector)pp));
    float   aa = sqrt(AREA((vector)pos));
    float   plimit = fw/lacun;
    float   scale = 1;
    int     octave = 0;
    float   fbm = 0;
    VOP_DAMPENFBM_FLOAT
    fbm = amp * fit(fbm, -1.25, 1.25, 0, 1);
    return fbm;
}

#define VOP_DAMPENFBM_VECTOR \
        while (scale > plimit && octave < maxoctaves) { \
        fbm += scale * vector(fit(vector(noise(pp)),{0.3,0.3,0.3},{0.7,0.7,0.7},{-1.0,-1.0,-1.0},{1.0,1.0,1.0})); \
        scale *= lacun; \
        pp /= lacun; \
        octave ++; \
        } \
        scale *= lacun; \
        fbm += scale * vector(fit(vector(noise(pp)),{0.3,0.3,0.3},{0.7,0.7,0.7},{-1.0,-1.0,-1.0},{1.0,1.0,1.0})); \
        while (scale > aa) { \
        fbm += scale * vector(fit(vector(noise(pp)),{0.3,0.3,0.3},{0.7,0.7,0.7},{-1.0,-1.0,-1.0},{1.0,1.0,1.0})); \
        pp *= lacun; \
        scale *= lacun; \
	    } \

vector
vop_dampenFbmVF(float pos, freq, offset; float lacun, amp; int maxoctaves)
{
    float  pp = pos * freq + offset;
    float  fw = vop_FilterWidth(pp);
    float  aa = vop_FilterWidth(pos);
    float  plimit = fw/lacun;
    float  scale = 1;
    int    octave = 0;
    vector fbm = 0;
    VOP_DAMPENFBM_VECTOR
    fbm = amp * vector(fit(fbm,{-1.25,-1.25,-1.25},{1.25,1.25,1.25},{0,0,0},{1,1,1}));
    return fbm;
}

vector
vop_dampenFbmVV(vector pos, freq, offset; float lacun, amp; int maxoctaves)
{
    vector pp = pos * freq + offset;
    float  fw = sqrt(AREA(pp));
    float  aa = sqrt(AREA(pos));
    float  plimit = fw/lacun;
    float  scale = 1;
    int    octave = 0;
    vector fbm = 0;
    VOP_DAMPENFBM_VECTOR
    fbm = amp * vector(fit(fbm,{-1.25,-1.25,-1.25},{1.25,1.25,1.25},{0,0,0},{1,1,1}));
    return fbm;
}

vector
vop_dampenFbmVP(vector4 pos, freq, offset; float lacun, amp; int maxoctaves)
{
    vector4	pp = pos * freq + offset;
    float fw = sqrt(AREA((vector)pp));
    float aa = sqrt(AREA((vector)pos));
    float plimit = fw/lacun;
    float scale = 1;
    int   octave = 0;
    vector	fbm = 0;
    VOP_DAMPENFBM_VECTOR
    fbm = amp * vector(fit(fbm,{-1.25,-1.25,-1.25},{1.25,1.25,1.25},{0,0,0},{1,1,1}));
    return fbm;
}

#define VOP_GENERICNOISE_FUNC(type, noisefunc, ZERO, SCALE) \
	    float	scale	= amp; \
	    int		i; \
	    nval = 0; \
	    for (i = 0; i < turb; i++, pp *= 2.0, scale *= rough) \
		nval += SCALE * scale * ((type(noisefunc(pp))) + ZERO); \
	    nval = (type(pow(nval, atten)));

#define VOP_PERLINNOISE_FUNC(type, ZERO, SCALE)	\
    VOP_GENERICNOISE_FUNC(type, noise, ZERO, SCALE)

#define VOP_SIMPLEXNOISE_FUNC(type)	\
    VOP_GENERICNOISE_FUNC(type, xnoise, -.5, .5)

#define VOP_PERLINCURLNOISE_FUNC(type) \
    VOP_GENERICNOISE_FUNC(type, curlnoise, 0, 1)

#define VOP_SIMPLEXCURLNOISE_FUNC(type) \
    VOP_GENERICNOISE_FUNC(type, curlxnoise, 0, 1)

#define VOP_SIMPLEXCURLNOISE2D_FUNC(type) \
    VOP_GENERICNOISE_FUNC(type, curlxnoise2d, 0, 1)

#define VOP_PERLINCURLNOISE2D_FUNC(type) \
    VOP_GENERICNOISE_FUNC(type, curlnoise2d, 0, 1)

float
vop_perlinNoiseVF(vector pos; int turb; float amp, rough, atten)
{
    vector pp = pos;
    float nval;
    VOP_PERLINNOISE_FUNC(float, 0, 1)
    return nval;
}

vector
vop_perlinNoiseVV(vector pos; int turb; float amp, rough, atten)
{
    vector pp = pos;
    vector nval;
    VOP_PERLINNOISE_FUNC(vector, 0, 1)
    return nval;
}

vector
vop_perlinNoiseVP(vector4 pos; int turb; float amp, rough, atten)
{
    vector4 pp = pos;
    vector nval;
    VOP_PERLINNOISE_FUNC(vector, 0, 1)
    return nval;
}

float
vop_correctperlinNoiseVF(vector pos; int turb; float amp, rough, atten)
{
    vector pp = pos;
    float nval;
    VOP_PERLINNOISE_FUNC(float, -.5, 0.5)
    return nval;
}

vector
vop_correctperlinNoiseVV(vector pos; int turb; float amp, rough, atten)
{
    vector pp = pos;
    vector nval;
    VOP_PERLINNOISE_FUNC(vector, -.5, 0.5)
    return nval;
}

vector
vop_correctperlinNoiseVP(vector4 pos; int turb; float amp, rough, atten)
{
    vector4 pp = pos;
    vector nval;
    VOP_PERLINNOISE_FUNC(vector, -.5, 0.5)
    return nval;
}

float
vop_simplexNoiseVF(vector pos; int turb; float amp, rough, atten)
{
    vector pp = pos;
    float nval;
    VOP_SIMPLEXNOISE_FUNC(float)
    return nval;
}

vector
vop_simplexNoiseVV(vector pos; int turb; float amp, rough, atten)
{
    vector pp = pos;
    vector nval;
    VOP_SIMPLEXNOISE_FUNC(vector)
    return nval;
}

vector
vop_simplexNoiseVP(vector4 pos; int turb; float amp, rough, atten)
{
    vector4 pp = pos;
    vector nval;
    VOP_SIMPLEXNOISE_FUNC(vector)
    return nval;
}

vector
vop_simplexCurlNoiseVV(vector pos; int turb; float amp, rough, atten)
{
    vector pp = pos;
    vector nval;
    VOP_SIMPLEXCURLNOISE_FUNC(vector)
    return nval;
}

vector
vop_simplexCurlNoise2DVV(vector pos; int turb; float amp, rough, atten)
{
    vector pp = pos;
    vector nval;
    VOP_SIMPLEXCURLNOISE2D_FUNC(vector)
    return nval;
}

vector
vop_simplexCurlNoiseVP(vector4 pos; int turb; float amp, rough, atten)
{
    vector4 pp = pos;
    vector nval;
    VOP_SIMPLEXCURLNOISE_FUNC(vector)
    return nval;
}

vector
vop_perlinCurlNoiseVV(vector pos; int turb; float amp, rough, atten)
{
    vector pp = pos;
    vector nval;
    VOP_PERLINCURLNOISE_FUNC(vector)
    return nval;
}

vector
vop_perlinCurlNoise2DVV(vector pos; int turb; float amp, rough, atten)
{
    vector pp = pos;
    vector nval;
    VOP_PERLINCURLNOISE2D_FUNC(vector)
    return nval;
}

vector
vop_perlinCurlNoiseVP(vector4 pos; int turb; float amp, rough, atten)
{
    vector4 pp = pos;
    vector nval;
    VOP_PERLINCURLNOISE_FUNC(vector)
    return nval;
}

#define VOP_CURLNOISE_SCALE_POTENTIAL(p, potential)		\
    if( geo != "" )						\
    {								\
	dist = volumesample(geo, 0, (vector)p);			\
	r = abs( clamp(dist/radius, -1, 1) );			\
	norm = normalize( volumegradient(geo, 0, (vector)p) );	\
	if( dist < 0 )						\
	{							\
	    norm = -norm;					\
	}							\
    }								\
    else							\
    {								\
	r = abs( clamp(distance/radius, -1, 1) );		\
    }								\
    r = ( 15.0*r - 10.0*r*r*r + 3*r*r*r*r*r ) / 8.0;		\
    potential = (r * potential) + ((1-r) * dot(norm, potential) * norm);

#define VOP_CURLNOISE_FUNC()	\
    float dist, r, d;					\
    vector norm = normalize(nml);			\
    VOP_CURLNOISE_SCALE_POTENTIAL(pos, noisevec);	\
    VOP_CURLNOISE_SCALE_POTENTIAL(xDiff, xDiffNoise);	\
    VOP_CURLNOISE_SCALE_POTENTIAL(yDiff, yDiffNoise);	\
    VOP_CURLNOISE_SCALE_POTENTIAL(zDiff, zDiffNoise);	\
    float dzdy = (yDiffNoise.z - noisevec.z);		\
    float dydz = (zDiffNoise.y - noisevec.y);		\
    float dxdz = (zDiffNoise.x - noisevec.x);		\
    float dzdx = (xDiffNoise.z - noisevec.z);		\
    float dydx = (xDiffNoise.y - noisevec.y);		\
    float dxdy = (yDiffNoise.x - noisevec.x);		\
    val.x = (dzdy - dydz)/h;				\
    val.y = (dxdz - dzdx)/h;				\
    val.z = (dydx - dxdy)/h;				\
    if( bounce )					\
    {							\
	vector vn;					\
	if( geo != "" )					\
	{						\
	    dist = volumesample(geo, 0, (vector)pos);	\
	    if( dist < 0 )				\
	    {						\
		norm = normalize( volumegradient(geo, 0, (vector)pos) ); \
		d = dot(norm, val);			\
		if( d < 0 )				\
		{					\
		    vn = d * norm;			\
		    val = val - 2*vn;			\
		}					\
	    }						\
	}						\
	else if( distance < 0 )				\
	{						\
	    d = dot(norm, val);				\
	    if( d < 0 )					\
	    {						\
		vn = d * nml;			\
		val = val - 2*vn;				\
	    }						\
	}						\
    }							
    
vector
vop_curlNoiseVV(vector pos, freq, offset, nml; 
		string type, geo;
		int turb, bounce;
	    	float amp, rough, atten, distance, radius, h)
{
    vector val = {0,0,0};
    
    if (type == "exact_pnoise")
    {
	return vop_perlinCurlNoiseVV(pos*freq-offset, turb, amp, rough*2, atten);
    }
    else if (type == "exact_xnoise")
    {
	return vop_simplexCurlNoiseVV(pos*freq-offset, turb, amp, rough*2, atten);
    }

    // Finite difference helpers
    vector xDiff = pos;		xDiff.x += h;
    vector yDiff = pos;		yDiff.y += h;
    vector zDiff = pos;		zDiff.z += h;

    // Noise vectors at and around pos
    vector noisevec, xDiffNoise, yDiffNoise, zDiffNoise;
    if( type == "onoise" )
    {
	noisevec = onoise(pos*freq - offset, turb, rough, atten) * amp;
	xDiffNoise = onoise(xDiff*freq - offset, turb, rough, atten) * amp;
	yDiffNoise = onoise(yDiff*freq - offset, turb, rough, atten) * amp;
	zDiffNoise = onoise(zDiff*freq - offset, turb, rough, atten) * amp;
    }
    else if( type == "snoise" )
    {
	noisevec = snoise(pos*freq - offset, turb, rough, atten) * amp;
	xDiffNoise = snoise(xDiff*freq - offset, turb, rough, atten) * amp;
	yDiffNoise = snoise(yDiff*freq - offset, turb, rough, atten) * amp;
	zDiffNoise = snoise(zDiff*freq - offset, turb, rough, atten) * amp;
    }
    else if( type == "anoise" )
    {
	noisevec = anoise(pos*freq - offset, turb, rough, atten) * amp;
	xDiffNoise = anoise(xDiff*freq - offset, turb, rough, atten) * amp;
	yDiffNoise = anoise(yDiff*freq - offset, turb, rough, atten) * amp;
	zDiffNoise = anoise(zDiff*freq - offset, turb, rough, atten) * amp;
    }
    else if (type == "xnoise" )
    {
	noisevec = vop_simplexNoiseVV(pos*freq - offset, turb, amp, rough, atten);
	xDiffNoise = vop_simplexNoiseVV(xDiff*freq - offset, turb, amp, rough, atten);
	yDiffNoise = vop_simplexNoiseVV(yDiff*freq - offset, turb, amp, rough, atten);
	zDiffNoise = vop_simplexNoiseVV(zDiff*freq - offset, turb, amp, rough, atten);
    }
    else
    {
	noisevec = vop_perlinNoiseVV(pos*freq - offset, turb, amp, rough, atten);
	xDiffNoise = vop_perlinNoiseVV(xDiff*freq - offset, turb, amp, rough, atten);
	yDiffNoise = vop_perlinNoiseVV(yDiff*freq - offset, turb, amp, rough, atten);
	zDiffNoise = vop_perlinNoiseVV(zDiff*freq - offset, turb, amp, rough, atten);
    }

    // Ramp function, partial derivatives and cross product
    VOP_CURLNOISE_FUNC()

    return val;
}

vector
vop_curlNoiseVP(vector4 pos, freq, offset; 
		vector nml; 
		string type; string geo;
		int turb, bounce;
		float amp, rough, atten, distance, radius, h)
{
    vector val = {0,0,0};

    if (type == "exact_pnoise")
    {
	return vop_perlinCurlNoiseVP(pos*freq-offset, turb, amp, rough*2, atten);
    }
    else if (type == "exact_xnoise")
    {
	return vop_simplexCurlNoiseVP(pos*freq-offset, turb, amp, rough*2, atten);
    }

    // Finite difference helpers
    vector4 xDiff = pos;	xDiff.x += h;
    vector4 yDiff = pos;	yDiff.y += h;
    vector4 zDiff = pos;	zDiff.z += h;

    vector noisevec, xDiffNoise, yDiffNoise, zDiffNoise;

    // Noise vectors at and around pos
    if (type == "xnoise")
    {
	noisevec = vop_simplexNoiseVP(pos*freq - offset, turb, amp, rough, atten);
	xDiffNoise = vop_simplexNoiseVP(xDiff*freq - offset, turb, amp, rough, atten);
	yDiffNoise = vop_simplexNoiseVP(yDiff*freq - offset, turb, amp, rough, atten);
	zDiffNoise = vop_simplexNoiseVP(zDiff*freq - offset, turb, amp, rough, atten);
    }
    else
    {
	noisevec = vop_perlinNoiseVP(pos*freq - offset, turb, amp, rough, atten);
	xDiffNoise = vop_perlinNoiseVP(xDiff*freq - offset, turb, amp, rough, atten);
	yDiffNoise = vop_perlinNoiseVP(yDiff*freq - offset, turb, amp, rough, atten);
	zDiffNoise = vop_perlinNoiseVP(zDiff*freq - offset, turb, amp, rough, atten);
    }

    // Ramp function, partial derivatives and cross product
    VOP_CURLNOISE_FUNC()

    return val;
}

vector
vop_curlNoise2DVV(vector pos, freq, offset;
		    string type, geo;
		    int turb;
		    float amp, rough, atten, distance, radius, h)
{
    vector val = {0,0,0};

    if (type == "exact_pnoise")
    {
	return vop_perlinCurlNoise2DVV(pos*freq-offset, turb, amp, rough*2, atten);
    }
    else if (type == "exact_xnoise")
    {
	return vop_simplexCurlNoise2DVV(pos*freq-offset, turb, amp, rough*2, atten);
    }

    // Finite difference helpers
    vector xDiff = pos;		xDiff.x += h;
    vector yDiff = pos;		yDiff.y += h;

    // Noise vectors at and around pos
    float noise, xDiffNoise, yDiffNoise;
    if( type == "onoise" )
    {
	noise = onoise(pos*freq - offset, turb, rough, atten) * amp;
	xDiffNoise = onoise(xDiff*freq - offset, turb, rough, atten) * amp;
	yDiffNoise = onoise(yDiff*freq - offset, turb, rough, atten) * amp;
    }
    else if( type == "snoise" )
    {
	noise = snoise(pos*freq - offset, turb, rough, atten) * amp;
	xDiffNoise = snoise(xDiff*freq - offset, turb, rough, atten) * amp;
	yDiffNoise = snoise(yDiff*freq - offset, turb, rough, atten) * amp;
    }
    else if( type == "anoise" )
    {
	noise = anoise(pos*freq - offset, turb, rough, atten) * amp;
	xDiffNoise = anoise(xDiff*freq - offset, turb, rough, atten) * amp;
	yDiffNoise = anoise(yDiff*freq - offset, turb, rough, atten) * amp;
    }
    else if (type == "xnoise")
    {
	noise = vop_simplexNoiseVF(pos*freq - offset, turb, amp, rough, atten);
	xDiffNoise = vop_simplexNoiseVF(xDiff*freq - offset, turb, amp, rough, atten);
	yDiffNoise = vop_simplexNoiseVF(yDiff*freq - offset, turb, amp, rough, atten);
    }
    else
    {
	noise = vop_perlinNoiseVF(pos*freq - offset, turb, amp, rough, atten);
	xDiffNoise = vop_perlinNoiseVF(xDiff*freq - offset, turb, amp, rough, atten);
	yDiffNoise = vop_perlinNoiseVF(yDiff*freq - offset, turb, amp, rough, atten);
    }

    // Scale the noise field for collision geometry.
    float r;
    if( geo != "" )
    {
	float dist = volumesample(geo, 0, pos);
	r = clamp(dist/radius, -1, 1);
    }
    else
    {
	r = clamp(distance/radius, -1, 1);
    }
    r = ( 15.0 * r - 10.0 * r*r*r + 3.0 * r*r*r*r*r ) / 8.0;
    noise = noise * r;
    xDiffNoise = xDiffNoise * r;
    yDiffNoise = yDiffNoise * r;

    // Take the curl of the noise field
    val.x = (yDiffNoise - noise) / h;
    val.y = (noise - xDiffNoise) / h;

    return val;
}

void
vop_displaceAlongNormal(vector pp, nn; float amount, scale, sshear, tshear;
			    int	 for_poly, obj_space, bump_only;
			    vector dP, dN)
{
#if defined(VOP_SHADING)
    if (obj_space)
    {
	float	nscale = length(nn);
	dN = normalize(ntransform("space:object", nn)) * nscale;
	dP = ptransform("space:object", pp);
    }
    else
    {
	dN = nn;
	dP = pp;
    }

    dP += (scale * amount)*dN;
    dP += (sshear * amount) * normalize(Du(dP));
    dP += (tshear * amount) * normalize(Dv(dP));

    if (obj_space)
    {
	dP = ptransform("space:object", "space:current", dP);
    }

    if (for_poly)
	 dN = computenormal(dP, nn, Ng);
    else dN = computenormal(dP);
#else
    dN = nn;
    dP = pp;
    dP += (scale * amount)*dN;
#endif

    if (bump_only)
    {
	dP =  pp;
    }
}

/***********************************************************************************/
/*  This is NOT working correctly. Will need to be reworked when soemone has time  */
/***********************************************************************************/
#if defined(VOP_SHADING)
void
vop_displaceAlongVector(vector pp, vec; string space; int mode;
                        float amount, scale, sshear, tshear;
			            vector dP, dN)
{
    dP = pp;
    dN = normalize(N);

    // transform into desired space and normalize vector
    vector ng = Ng;
    vector oP = pp;
    vector oN = vec;
    if (space != "space:current" ) {

        oP = ptransform(space, pp);
        ng = ntransform(space, Ng);

        if ( mode == 0 )                // normal mapping
            oN = amount * normalize(ntransform(space, 2*vec-1));
        else if ( mode == 1 )           // displace along vector use magnitude as amount
            oN = ntransform(space, vec);
        else if ( mode == 2 )           // displace along normal
            oN = amount * normalize(ntransform(space, N));
    }
    else {
        if ( mode == 0 )                // normal mapping
            oN = amount * normalize(2*vec-1);
        else if ( mode == 1 )           // displace along vector use magnitude as amount
            oN = vec;
        else if ( mode == 2 )           // displace along normal
            oN = amount * normalize(N);
    }
    
    // do displacement
    oP += scale * oN;
	oN = computenormal(oP, oN, ng);


    // return to original space
    if (space != "space:current" ) {
        dP = ptransform(space,"space:current", oP);
        dN = normalize(ntransform(space,"space:current", oN));
    }
    else {
        dP = oP;
        dN = normalize(oN);
    }
}
#endif

vector
vop_setcomp(vector in; float fval; int part)
{
    vector out = in;
    setcomp(out, fval, part);
    return out;
}

matrix
vop_setmatcomp(matrix in; float fval; int row, col)
{
    matrix out = in;
    setcomp(out, fval, row, col);
    return out;
}

#define VOP_TRANSLATE_FN(M_TYPE, V_TYPE) \
M_TYPE \
vop_translate(M_TYPE in; V_TYPE t) \
{ \
    M_TYPE out = in; \
    translate(out, t); \
    return out; \
} 

VOP_TRANSLATE_FN(matrix, vector)
VOP_TRANSLATE_FN(matrix, vector4)
#undef VOP_TRANSLATE_FN

#define VOP_ROTATE_FN(M_TYPE) \
M_TYPE \
vop_rotate(M_TYPE in; float angle; vector axis) \
{ \
    M_TYPE out = in; \
    rotate(out, angle, axis); \
    return out; \
}

VOP_ROTATE_FN(matrix)
VOP_ROTATE_FN(matrix3)
#undef VOP_ROTATE_FN


#define VOP_SCALE_FN(M_TYPE) \
M_TYPE \
vop_scale(M_TYPE in; vector s) \
{ \
    M_TYPE out = in; \
    scale(out, s); \
    return out; \
}

VOP_SCALE_FN(matrix)
VOP_SCALE_FN(matrix3)
#undef VOP_SCALE_FN

vector
vop_frompolar(float u, v; float radius)
{
    float	sv = sin(v);
    return set(sv*cos(u), sv*sin(u), cos(v))*radius;
}

vector
vop_topolarXYZ(float x, y, z)
{
    float r = sqrt(x*x+y*y+z*z);
    return set(
	    atan(y, x) % (2* M_PI),
	    acos(z/r),
	    r);
}

vector
vop_topolar(vector v)
{
    return vop_topolarXYZ(v.x, v.y, v.z);
}

#if defined(VOP_SHADING)
void
vop_computeTangents(string tstyle;
		    vector nn, uv;
		    vector in_utan, in_vtan;
		    vector out_utan, out_vtan)
{
    if (tstyle == "geo")
    {
	out_utan = normalize(dPds);
	out_vtan = normalize(dPdt);
    }
    else if (tstyle == "world")
    {
	out_vtan = cross(set(0, 0, 1), nn);
	out_vtan = length2(out_vtan) < 1e-6 ? set(1, 0, 0) : normalize(out_vtan);
	out_utan = normalize(cross(nn, out_vtan));
    }
    else if (tstyle == "object")
    {
	out_vtan = cross(ow_vspace(set(0, 0, 1)), nn);
	out_vtan = length2(out_vtan) < 1e-6 ? set(1, 0, 0) : normalize(out_vtan);
	out_utan = normalize(cross(nn, out_vtan));
    }
    else if (tstyle == "uv")
    {
	out_utan = normalize((dPds * Dv(uv.y)) - (dPdt * Du(uv.y)));
	out_vtan = normalize(cross(nn, out_utan));
	out_utan = normalize(cross(nn, out_vtan));
    }
    else
    {
	out_utan = in_utan;
	out_vtan = in_vtan;
    }
}

void
vop_computeTangentsOld(vector out_tanu, out_tanv;
		       vector nn, uv;
		       vector in_tanu, in_tanv;
		       int tstyle)
{
    string	s_tstyle;
    if (tstyle == 0)
	s_tstyle = "world";
    else if (tstyle == 1)
	s_tstyle = "object";
    else if (tstyle == 2)
    s_tstyle = "uv";
    else
	s_tstyle = "inputs";
    vop_computeTangents(s_tstyle,
	    nn, uv, in_tanu, in_tanv, out_tanu, out_tanv);
}

vector
vop_anisotropic_eval(vector ll; vector nn; vector V;
		vector uv; vector in_tanu; vector in_tanv;
		float urough, vrough; int model; int tstyle)
{
    vector    H;	// Half angle vector
    vector    clr;	// Color
    vector    tanU, tanV;
    float     rz, cos_r, cos_i; // Reflected and incident angles

    float     nml_term;
    float     uval, vval, nval;
    float     exponent;

    cos_r = vop_dot(nn, V);
    clr = 0;
    if (cos_r > 0.0)
    {
	vop_computeTangentsOld(tanU, tanV, nn, uv, in_tanu, in_tanv, tstyle);

	cos_i = vop_dot(ll, nn);
	if (cos_i > 0.0)
	{
	    H = normalize(V + ll);
	    uval = vop_dot(tanU, H);
	    vval = vop_dot(tanV, H);
	    nval = vop_dot(nn, H);

	    rz = 0;
	    if (nval > 0)
	    {
		if (model == 0)
		{
		    // Ward model
		    nml_term = 4.0 * M_PI * urough*vrough;
		    uval /= urough;
		    vval /= vrough;
		    rz = cos_i*exp(-2.*(uval*uval + vval*vval) /
			    (1.0 + nval)); 
		    rz /= nml_term * sqrt(cos_i*cos_r); 
		}
		else
		{
		    // Ashikhmin model
		    exponent = uval*uval/urough + vval*vval/vrough;
		    exponent /= 1.0 - nval*nval;

		    rz = pow(nval, exponent) / (4.0 * dot(V, H));
		}
	    }
	    clr = rz;
	}
    }
    return clr;
}
#endif

#if defined(VOP_SURFACE) || defined(VOP_DISPLACE) || defined(VOP_FOG)
vector
vop_anisotropic(vector nn; vector V; vector uv; float urough, vrough;
		int model, tstyle)
{
    vector	ll;              // Normalized light vector
    vector	lclr;            // Light color
    vector	tanu = 0;
    vector	tanv = 0;

    lclr = 0;
    tanu = 0;
    tanv = 0;

    // Loop through all lights which contribute specular
    illuminance (P, nn, M_PI/2, LIGHT_SPECULAR, "lightexport", "")
    {
	shadow(Cl);
	ll = normalize(L);
	lclr = vop_anisotropic_eval(
		ll, nn, V, uv, tanu, tanv, urough, vrough, model, tstyle);
	lclr *= Cl;
    }
    return lclr;
}

vector
vop_sheen(vector nn, ii; float eta, rough; int facefwd)
{
    vector	R, T;
    float	Kr, Kt, sheen;
	vector	nf = nn;
	if (facefwd) nf = normalize(vop_frontface(nn, ii));
    vector	illum = 0;

    fresnel(ii, nn, eta, Kr, Kt, R, T);
    Kr = smooth(0.0, 0.5, Kr);
    illuminance (P, nn, M_PI/2, LIGHT_SPECULAR, "lightexport", "")
    {
	vector nL = normalize(L);
	shadow(Cl);
	sheen = specularBRDF(nL, nf, -ii, rough);
	illum += Cl * vop_dot(nL, nf) * (sheen + 0.2);
    }
    return Kr * illum;
}
#endif

#if defined(VOP_SHADING)
bsdf
vop_sheen_bsdf(vector nn, ii; float eta, rough; int facefwd)
{
    float	Kr, Kt;
    vector	R, T;
	vector	nf = nn;
	if (facefwd) nf = normalize(vop_frontface(nn, ii));
	bsdf    f;

    fresnel(ii, nn, eta, Kr, Kt, R, T);
    Kr = smooth(0.0, 0.5, Kr);
    f = Kr * ((bsdf(diffuse(nf)) * 0.2) +
	       bsdf(diffuse(nf)) * matchvex_specular(nf, 1.0/rough));
    return f;
}

float
vop_specular_eval(string lmodel; vector ll, nf, ii, uv, tanu, tanv;
	     float urough, vrough, sharp; int tstyle)
{
    float	seval = 0;

    if (lmodel == "phong")
    {
	seval = phongBRDF(ll, nf, -ii, 1.0/urough);
    }
    else if (lmodel == "blinn")
    {
	seval = blinnBRDF(ll, nf, -ii, urough);
    }
    else if (lmodel == "glossy")
    {
	float	w = sharp/2.0;		// fit sharpness between 0 & 0.5

	seval = specularBRDF(ll, nf, -ii, urough); 
	if (w > 0)
	    seval = smooth(w, 1-w, seval); 
    }
    else if (lmodel == "anisotropic")
    {
	seval = (float)vop_anisotropic_eval(ll, nf, -ii,
		uv, tanu, tanv, urough, vrough, 0, tstyle);
    }
    else if (lmodel == "spec" ||
	     lmodel == "specular")
    {
	seval = specularBRDF(ll, nf, -ii, urough); 
    }
    return seval;
}
#endif

#if defined(VOP_SURFACE) || defined(VOP_DISPLACE) || defined(VOP_FOG)
vector
vop_specular(string lmodel; vector nf, ii, uv;
	     float urough, vrough, sharp; int tstyle)
{
    vector	clr, ll;
    vector	tanu, tanv;

    tanu = 0;
    tanv = 0;

    clr = 0;
    illuminance(P, nf, PI/2, LIGHT_SPECULAR, "lightexport", "")
    {
	shadow(Cl);
	ll = normalize(L);
	clr += Cl * vop_specular_eval(lmodel, ll, nf, ii, uv, tanu, tanv,
		urough, vrough, sharp, tstyle);
    }

    return clr;
}
#endif

#if defined(VOP_SHADING)
bsdf
vop_specular_bsdf(string lmodel; vector nf, ii, uv, tanu, tanv;
	     float urough, vrough, sharp; int tstyle)
{
    bsdf	f;

    if (lmodel == "phong")
    {
	f = bsdf(phong(nf, 1.0/urough));
    }
    else if (lmodel == "blinn")
    {
	f = matchvex_blinn(nf, 1.0/urough);
    }
    else if (lmodel == "glossy")
    {
	f = matchvex_specular(nf, 1.0/urough);
    }
    else if (lmodel == "anisotropic")
    {
	bsdf	aniso_bsdf;
	vector	tanU;
	vector	tanV;

	vop_computeTangentsOld(tanU, tanV, nf, uv, tanu, tanv, tstyle);
	aniso_bsdf = ashikhmin(nf, 2.0/(urough*urough), 2.0/(vrough*vrough),
			       normalize(tanU), normalize(tanV));
	f = (1.0 / (2.0 * M_PI * urough * vrough)) * aniso_bsdf;
    }
    else
    {
	f = matchvex_specular(nf, 1.0/urough);
    }
    return f;
}

float
vop_diffuse_eval(string dmodel; vector ll, nf, ii; float rough)
{
    float	deval = 0;
    if (dmodel == "isotropic")
    {
	deval = 1;
    }
    else if (dmodel == "oren")
    {
	deval = diffuseBRDF(ll, nf, -ii, rough);
    }
    else if (dmodel == "diffuse")
    {
	deval = diffuseBRDF(ll, nf);
    }
    return deval;
}

bsdf
vop_diffuse_bsdf(string dmodel; vector nf, ii; float rough)
{
    bsdf	f;

    if (dmodel == "isotropic")
    {
	f = isotropic();
    }
    else if (dmodel == "oren")
    {
	f = bsdf(diffuse(nf, rough));
    }
    else
    {
	f = bsdf(diffuse(nf));
    }
    return f;
}
#endif

#if defined(VOP_SURFACE) || defined(VOP_DISPLACE) || defined(VOP_FOG)
vector
vop_lighting(string lmodel; 
	     vector nf, ii, uv, amb, diff, spec; 
	     float urough, vrough; int tstyle)
{
    vector	clr;

    if (lmodel == "constant")
    {
	clr = diff;
    }
    else if (lmodel == "headlight")
    {
	clr = diff * clamp(vop_dot(nf, -ii), 0, 1);
    }
    else
    {
	float	angle;

	clr = amb * ambient();

	if (lmodel == "isotropic")
	    angle = PI;
	else
	    angle = PI/2;

	illuminance(P, nf, angle, LIGHT_DIFFSPEC, "lightexport", "")
	{
	    vector ll = normalize(L);
	    vector lclr = 0;

	    if (diff != {0,0,0})
	    {
		float	deval;
		if (lmodel == "isotropic" ||
		    lmodel == "oren")
		    deval = vop_diffuse_eval(lmodel, ll, nf, ii, urough);
		else
		    deval = vop_diffuse_eval("diffuse", ll, nf, ii, urough);

		lclr += deval * diff;
	    }
	    if (lmodel != "isotropic" &&
		lmodel != "oren" &&
		spec != {0,0,0})
	    {
		vector	tanu = 0;
		vector	tanv = 0;

		float seval = vop_specular_eval(lmodel, ll, nf, ii,
			uv, tanu, tanv, urough, vrough, 0, tstyle);

		lclr += seval * spec;
	    }

	    if (lclr != 0)
	    {
		shadow(Cl);
		clr += lclr * Cl;
	    }
	}
    }
    return clr;
}
#endif

#if defined(VOP_SHADING)
bsdf
vop_lighting_bsdf(string lmodel; 
	     vector nf, ii, uv, amb, diff, spec; 
	     float urough, vrough; int tstyle)
{
    bsdf	f;
    vector	tanu = 0;
    vector	tanv = 0;

    if (lmodel == "constant")
    {
	f = bsdf();
    }
    else if (lmodel == "lambert" ||
	     lmodel == "headlight" ||
	     lmodel == "oren" ||
	     lmodel == "isotropic")
    {
	f = diff * vop_diffuse_bsdf(lmodel, nf, ii, urough);
    }
    else
    {
	tanu = 0;
	tanv = 0;
	f = diff * bsdf(diffuse(nf));
	f += spec * vop_specular_bsdf(lmodel, nf, ii,
		uv, tanu, tanv, urough, vrough, 0, tstyle);
    }
    return f;
}
#endif

#if defined(VOP_SURFACE) || defined(VOP_DISPLACE) || defined(VOP_FOG)

#include "singlescatter.h"
#include "pcscatter.h"

vector
vop_hairspec(vector nn, V, T; float exp;)
{
    // Specular illumination model from:
    // James T. Kajiya and Timothy L.  Kay (1989) "Rendering Fur with Three 
    // Dimensional Textures", Computer Graphics 23,3, 271-280  

    float NdotI = vop_dot(nn, V);
    vector illum = 0;

    illuminance (P, nn, M_PI/2, LIGHT_SPECULAR, "lightexport", "") 
    {
	vector nL = normalize(L);
	float NdotL = vop_dot(nn, nL);

	float Kajiya = cos(abs(acos(dot(T,nL)) - acos(dot(-T,V))));

	shadow(Cl);
	illum += Cl * NdotL * NdotI * pow(Kajiya, 1.0/exp); 
    }
    return illum;
}

vector
vop_simpleSSS(vector nn, ii; float  eta, depth; int facefwd)
{
    // from Matt Pharr SIGGRAPH 2001 course notes
    // Based on phenomenological information about skin reflectance 
    // from Hanrahan and Krueger, proceedings of Siggraph 1993: 
    // "Reflection from layered surfaces due to subsurface scattering", 
    // the ratio of the indices of refraction of the incoming ray 
    // to the index of refraction for the transmitted ray
    // the overall depth of the surface layer to the subsurface layer
    // Loops over light sources and computes the reflected color.
    //
    // This performs a single iteration of sub-surface scattering.
    //
    vector	illum = 0;
    float	Kr, fKr = 1;
    float	Kt = 1;
    float	fKt = 1;
    vector	rayin = efresnel(ii, nn, eta, Kr, Kt);
	vector	nf = nn;
	if (facefwd) nf = normalize(vop_frontface(nn, ii));

    illuminance (P, nn, M_PI/2, LIGHT_DIFFUSE, "lightexport", "")
    {
	vector nL = normalize(L);
	vector rayout = efresnel(-nL, nn, eta, fKr, fKt);

	shadow(Cl);
	illum += Cl * vop_dot(nL, nf) * Kt * fKt * 
			    (singleScatter(rayin, rayout, nf, .8, .8, depth) +
			     singleScatter(rayin, rayout, nf, .3, .5, depth) +
			     singleScatter(rayin, rayout, nf, 0., .4, depth));
    }
    return illum;
}
#endif

#if defined(VOP_SHADING)
bsdf
vop_SSS_bsdf(vector nn, ii; float  eta, depth; int facefwd)
{
    float	Kr = 1;
    float	Kt = 1;
	vector  R, T;
	vector	nf = nn;
	if (facefwd) nf = normalize(vop_frontface(nn, ii));
	bsdf    f;

    fresnel(ii, nn, eta, Kr, Kt, R, T);
	Kr = smooth(0, 0.5, Kr);
	Kt = 1 - Kr;

	// this is incorrect, but the only current approximation. 23-07-2007
	f = Kt * bsdf(diffuse(nf));
    return f;
}
#endif

#if defined(VOP_SURFACE) || defined(VOP_DISPLACE) || defined(VOP_FOG)
vector
vop_multiSSS(vector Pin, Nin, Rd; float sd, bounce;
			 string pcmap; int nfp, t_rgb;)
{
   vector Xo = wo_space(Pin);
   vector No = normalize(wo_nspace(Nin));
   vector mapP, mapN, ssm;
   int xxx;
   string ch_ssm = "ssM";
   int handle = pcopen(pcmap, "P", Xo, "N", No, 1e37, nfp);
   while (pcunshaded(handle, ch_ssm)) {
      pcimport(handle, "P", mapP);
      pcimport(handle, "N", mapN);
      ssm = vop_ssIntegMulti(pcmap, Rd, sd, bounce,t_rgb, mapP, mapN);
      xxx = pcexport(handle, ch_ssm, ssm);
   }
   vector bssrdf = vector(pcfilter(handle, ch_ssm));
   pcclose(handle);
   return bssrdf;
}

vector 
vop_singleSSS(vector Pin, Nin, Iin, Rd; float sd, g, eta, tbias;
                string pcmap; int nfp, samples, t_rgb;)
{
   vector bssrdf;
   if(pcmap != "")
   {
      vector Xo = wo_space(Pin);
      vector No = normalize(wo_nspace(Nin));
      vector eye = Pin-Iin;
      vector mapP, mapN, ssm;
      int xxx;
      int handle = pcopen(pcmap, "P", Xo, "N", No, 1e37, nfp);
      string ch_ssm = "ssS";
      while (pcunshaded(handle, ch_ssm)) {
         pcimport(handle, "P", mapP);
         pcimport(handle, "N", mapN);
         ssm = vop_ssIntegSingle (Rd,sd, g,eta,samples,tbias,t_rgb,
                              ow_space(mapP),ow_nspace(mapN),mapP-eye);
         xxx = pcexport(handle, ch_ssm, ssm);
      }
      bssrdf = vector(pcfilter(handle, ch_ssm));
      pcclose(handle);
   }
   else
   {
      bssrdf = vop_ssIntegSingle(Rd,sd,g,eta,samples,tbias,t_rgb,Pin,Nin,Iin);
   }
   return bssrdf;
}
#endif

vector
vop_trace_environment(string envmap, envobj; vector raydir, bg;)
{
    vector val = 0;
    if (envmap == "")
    {
	val = bg;
    }
    else
    {
	vector	dir = vtransform("space:current", envobj, raydir);
	val = environment(envmap, dir, dir, dir, dir, "lerp", 1);
    }
    return val;
}

#if defined(VOP_SHADING)
vector
vop_trace(vector Nf, nN, PP, dir, atten, bg; string envmap;
	float bias, angle, thresh, jitter, density; string style;
	int samples; float nhit; string envobj;)
{
    nhit = 0;
    vector sum = 0;
    vector raydir = 0;
    float  atten_interp, dist;
    vector hitCf, hitOf, hitP;

    gather( PP, dir,
	    "Cf", hitCf,
	    "Of", hitOf,
	    "P",  hitP,
	    "ray:direction", raydir,
	    "samples", samples,
	    "bias", bias,
	    "angle", angle,
	    "raystyle", style,
	    "rayweight", thresh,
	    "samplebase", jitter,
	    "variancevar", "Cf")
    {
        // Attenuation
        atten_interp = 1;
        if ((dot(Nf, nN) > 0) && (dot(Nf, raydir) < 0) && (density > 0))
        {
            dist = length(hitP-PP);
            atten_interp = exp(-dist*density);
        }
        sum += lerp(atten, hitCf, atten_interp);

        if (vop_maxcomp(hitOf) < 1)
        {
            sum += (set(1,1,1) - hitOf) *
                    vop_trace_environment(envmap, envobj, raydir, bg);
        }

        nhit += 1;
    }
    else
    {
        // Resolve background color
        sum += vop_trace_environment(envmap, envobj, raydir, bg);
        nhit += 1;
    }
    return sum;
}
#endif

float
vop_weave(float ss, tt, sfreq, tfreq, width)
{
	float weave = 0;
	if ( sfreq > 0 && tfreq > 0 && width > 0)
	{
		// setup U V locations
		float vv = tt * 0.5;
		float uu = (ss+vv) * sfreq;
		vv = (ss-vv) * tfreq;
		if (vop_is_even(uu) ) vv += 0.5;
		
		// setup the pulse width
		float warppulse = sqrt(abs(sin(uu * M_PI)));
		warppulse = smooth(1-width, 1, warppulse);
		
		float weftpulse = sqrt(abs(sin(vv * M_PI)));
		weftpulse = smooth(1-width, 1, weftpulse);
		
		/* up and down in warp */
		float weft = lerp(0, 1, warppulse);
		weft = lerp(weft, 0, 1-warppulse);
		/* up and down in weft on warp */
		weft = lerp(0, weft, weftpulse);
		weft = lerp(weft, 0, 1-weftpulse);

		/* up and down in weft */
		float warp = lerp(0, 1, weftpulse);
		warp = lerp(warp, 0, 1-weftpulse);
		/* up and down in warp on weft */
		warp = lerp(0, warp, warppulse);
		warp = lerp(warp, 0, 1-warppulse);

		weave = ((vop_is_odd(uu) && vop_is_even(vv)) || 
						(vop_is_even(uu) && vop_is_odd(vv))) ?  
										weft : warp;
	}
	return weave;
}

vector
vop_toUnitNormal(vector Ni;)
{
    return Ni * 0.5 + 0.5;
}

vector
vop_fromUnitNormal(vector Ni;)
{
    return Ni * 2.0 - 1.0;
}

vector
vop_tangentNormal(vector Ni;
                  vector nn;
                  vector utan;
                  vector vtan;
                  int onspace;
                  int flipX;
                  int flipY;
                  float heightScale;)
{
    vector Nn   = normalize(nn);
    vector Nin  = lerp(Nn, Ni, heightScale);

    vector No;
    No.x = dot(Nin, utan);
    No.y = dot(Nin, vtan);
    No.z = dot(Nin, nn);
    No = normalize(No);

    if (flipX) { No.x = -No.x; }
    if (flipY) { No.y = -No.y; }

    if (onspace == 0) { No = vop_toUnitNormal(No); }

    return No;
}

vector
vop_tangentNormalRemap(vector Ni;
                       vector nn;
                       vector utan;
                       vector vtan;
                       int inspace;
                       int flipX;
                       int flipY;
                       float heightScale;)
{
    vector Nn   = normalize(nn);
    vector Vs	= normalize(utan);
    vector Vt	= normalize(vtan);
    vector Nin  = Ni;

    if (inspace == 0) { Nin = vop_fromUnitNormal(Nin); }

    if (flipX) { Nin.x = -Nin.x; }
    if (flipY) { Nin.y = -Nin.y; }

    matrix M;
    setcomp(M, Vs.x, 0, 0);
    setcomp(M, Vs.y, 0, 1);
    setcomp(M, Vs.z, 0, 2);
    setcomp(M, 0,    0, 3);
    setcomp(M, Vt.x, 1, 0);
    setcomp(M, Vt.y, 1, 1);
    setcomp(M, Vt.z, 1, 2);
    setcomp(M, 0,    1, 3);
    setcomp(M, Nn.x, 2, 0);
    setcomp(M, Nn.y, 2, 1);
    setcomp(M, Nn.z, 2, 2);
    setcomp(M, 0,    2, 3);
    setcomp(M, 0,    3, 0);
    setcomp(M, 0,    3, 1);
    setcomp(M, 0,    3, 2);
    setcomp(M, 1,    3, 3);

    vector No = normalize(ntransform(Nin, M));

    No = lerp(Nn, No, heightScale);

    return No;
}

vector
vop_bumpToNormalMap(string map;
                    int onspace;
                    int flipX;
                    int flipY;
                    float heightScale;
                    vector uv;)
{
    int xres = 0;
    teximport(map, "texture:xres", xres);
    int yres = 0;
    teximport(map, "texture:yres", yres);
    float dx = 1.0 / max(xres, yres);

    vector v00 = texture(map, uv.x-dx, uv.y-dx); // top left
    vector v01 = texture(map, uv.x-dx, uv.y   ); // left
    vector v02 = texture(map, uv.x-dx, uv.y+dx); // bottom left
    vector v10 = texture(map, uv.x,    uv.y-dx); // top
    vector v12 = texture(map, uv.x,    uv.y+dx); // bottom
    vector v20 = texture(map, uv.x+dx, uv.y-dx); // top right
    vector v21 = texture(map, uv.x+dx, uv.y   ); // right
    vector v22 = texture(map, uv.x+dx, uv.y+dx); // bottom right

    float f00 = luminance(v00);
    float f01 = luminance(v01);
    float f02 = luminance(v02);
    float f10 = luminance(v10);
    float f12 = luminance(v12);
    float f20 = luminance(v20);
    float f21 = luminance(v21);
    float f22 = luminance(v22);

    vector No;
    No.x = f20 + 2*f21 + f22 -f00 - 2*f01 - f02;
    No.y = f02 + 2*f12 + f22 -f00 - 2*f10 - f20;
    No.z = 1.0 / heightScale;
    if (flipX) { No.x = -No.x; }
    if (flipY) { No.y = -No.y; }
    No = normalize(No);

    if (onspace == 0) { No = vop_toUnitNormal(No); }

    return No;
}

#if defined(VOP_SHADING)
void
vop_curvature(vector p;
              vector n;
              int mode;
              int space;
              int smooth;
              float tolerance;
              float convexscale;
              float convexbias;
              float concavescale;
              float concavebias;
              float biasmap;
              float Ko;)
{
    vector xso = {1, 0, 0};
    vector xsw = ptransform("space:object", "space:world", xso);
    float xs = length(xsw);
    float xa = area(P);

    vector dPdu, dPdv;
    getderiv(p, "P", 0, s, t, dPdu, dPdv, "smooth", smooth);
    dPdu /= xa;
    dPdv /= xa;

    vector dNdu, dNdv;
    getderiv(n, "N", 0, s, t, dNdu, dNdv, "smooth", smooth);
    dNdu /= xa;
    dNdv /= xa;

    float a00 = dot(dPdu, dPdu);
    float a01 = dot(dPdu, dPdv);
    float a11 = dot(dPdv, dPdv);

    float b00 = -dot(dNdu, dPdu);
    float b01 = -(dot(dNdu, dPdv) + dot(dNdv, dPdu)) / 2;
    float b11 = -dot(dNdv, dPdv);

    float a = a00*a11 - a01*a01;
    float b = b00*b11 - b01*b01;
    float aa = max(abs(a), tolerance) * sign(a);

    float cn00 = a11 / aa;
    float cn01 = -a01 / aa;
    float cn11 = a00 / aa;

    float K = b / aa;
    float H = ((cn00*b00 + cn11*b11) / 2) + (cn01*b01);

    float Q = H*H - 4*K;
    float QQ = max(Q, 0);
    float SQ = sqrt(QQ);

    float p0 = (H + SQ) / 2;
    float p1 = (H - SQ) / 2;

    Ko = 0;

    if (mode == 0)
    {
        Ko = K;
    }
    else if (mode == 1)
    {
        Ko = -H;
    }

    Ko *= xs;

    float biasr = clamp(biasmap, 0, 1) - 0.5;
    if (Ko >  0) { Ko *= convexscale;  Ko =  vop_bias( Ko, clamp(convexbias  + biasr, 0, 1)); }
    if (Ko <= 0) { Ko *= concavescale; Ko = -vop_bias(-Ko, clamp(concavebias + biasr, 0, 1)); }

    if (space == 0)
    {
        Ko = fit(Ko, -1, 1, 0, 1);
    }
}
#endif

#endif
// :vimvfl:
