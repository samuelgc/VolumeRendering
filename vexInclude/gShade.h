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
 * NAME:	gShade.h ( Gallery Library, C++)
 *
 * COMMENTS:	Some useful "macros" for shading functions
 */

#ifndef __gShade__
#define __gShade__

#include "gBump.h"

//
// Compute shading space for 3D textures
//
#define gTransformRest(pp, restname, rest, space)	{ \
		 if (isbound(restname))	pp = rest; \
	    else if (space == "world")	pp = P; \
	    else if (space == "NDC")	pp = toNDC(P); \
	    else if (space == "object") pp = wo_space(P); \
	    else			pp = wt_space(P); }

//
// Compute UV texture coordinates for 2D textures
//
#define	gTextureCoord2D(ss, tt, name, uv)	 \
	    if (isbound(name))	{ ss = uv.x; tt = uv.y; } \
	    else		{ ss = s; tt = t; }

#define	gTextureCoord3D(st, name, uv)	 \
	    if (isbound(name))	{ st = uv; } \
	    else		{ st = set(s, t, .5); }

#define G_AMB_SHADEPARMS(prefix) \
    vector	prefix##ambient=1;

#define G_DIFF_SHADEPARMS(prefix) \
    vector	prefix##diffuse=.8;

#define G_DIFF_COLORPARMS(prefix, r, g, b) \
    vector	prefix##diffuse = CONST_COLOR(r, g, b);

#define G_SPEC_SHADEPARMS(prefix) \
    vector	prefix##specular =.8; \
    float	prefix##srough   = 0.05;

#define G_SPEC_COLORPARMS(prefix, r, g, b, rough) \
    vector	prefix##specular = CONST_COLOR(r, g, b); \
    float	prefix##srough   = rough;

#define G_REFL_SHADEPARMS(prefix) \
    vector	prefix##reflect = 0;

#define G_REFL_COLORPARMS(prefix, r, g, b) \
    vector	prefix##reflect = CONST_COLOR(r, g, b);

#define G_LMOD_SHADEPARMS(prefix) \
    string	prefix##lmodel   = "specular";

#define G_LMOD_COLORPARMS(prefix, model) \
    string	prefix##lmodel   = model;

#define G_BIAS_SHADERPARM()	\
    float	refl_bias = 0.01;

vector
gComputeSpecular(vector nn, vv; float rough; string lmodel)
{
    vector	spec = 0;

	 if (lmodel == "specular") spec = specular(nn, vv, rough);
    else if (lmodel == "blinn")	   spec = blinn(nn, vv, rough);
    else if (lmodel == "phong")	   spec = phong(nn, vv, 1.0/rough);
    return spec;
}

bsdf
gComputeSpecularBSDF(float rough; string lmodel)
{
    bsdf	spec = bsdf();

	 if (lmodel == "specular") spec = matchvex_specular(1.0/rough);
    else if (lmodel == "blinn")	   spec = matchvex_blinn(1.0/rough);
    else if (lmodel == "phong")	   spec = phong(1.0/rough);
    return spec;
}


vector
gComputeLight(vector nn, vv; vector amb, diff, spec, refl;
		float rough, refl_bias;
		string lmodel) 
{
    vector	clr;
    float	maxr;
    if (lmodel == "constant") clr = diff;
    else
    {
	clr  = diff*(amb*ambient() + diffuse(nn));
	clr += spec*gComputeSpecular(nn, vv, rough, lmodel);
    }
    maxr = (float)abs(refl);
    if (maxr > 0.0001)
	clr += refl * reflectlight(P, nn, -vv, refl_bias, maxr);
    return clr;
}

vector
gBlendTwoMaterials(vector nn, vv;
	    vector amb1, diff1, spec1, refl1; float rough1; string lmodel1;
	    vector amb2, diff2, spec2, refl2; float rough2; string lmodel2;
	    float refl_bias; float blend_amount)
{
    vector	clr, spec, aspec, bspec;
    float	maxr;
    if (lmodel1 == lmodel2)
    {
	if (lmodel1 == "constant")
	    clr = lerp(diff1, diff2, blend_amount);
	else
	{
	    clr = lerp(amb1, amb2, blend_amount)*ambient();
	    clr = lerp(diff1, diff2, blend_amount)*(clr + diffuse(nn));
	    if (lmodel1 != "lambert")
	    {
		if (rough1 == rough2)
		{
		    spec = gComputeSpecular(nn, vv, rough1, lmodel1);
		    spec *= lerp(spec1, spec2, blend_amount);
		}
		else
		{
		    aspec = gComputeSpecular(nn, vv, rough1, lmodel1);
		    bspec = gComputeSpecular(nn, vv, rough2, lmodel2);
		    spec = lerp(aspec*spec1, bspec*spec2, blend_amount);
		}
		clr += spec;
	    }
	}
    }
    else
    {
	// Here we're blending two separate lighting models.  So we sum up the
	// light without reflections first.
	aspec = gComputeLight(nn, vv, amb1, diff1, spec1, {0,0,0}, rough1,
				refl_bias, lmodel1);
	bspec = gComputeLight(nn, vv, amb2, diff2, spec2, {0,0,0}, rough2,
				refl_bias, lmodel2);
	clr = lerp(aspec, bspec, blend_amount);
    }
    spec = lerp(refl1, refl2, blend_amount);
    maxr = (float)abs(spec);
    if (maxr > 0.0001)
	clr += spec * reflectlight(P, nn, -vv, refl_bias, maxr);
    return clr;
}

bsdf
gComputeBSDF(vector nn, vv; vector amb, diff, spec, refl;
		float rough, refl_bias;
		string lmodel) 
{
    bsdf	fval;
    float	maxr;

    if (lmodel == "constant")
	fval = bsdf();
    else
	fval = diff * diffuse() + spec * gComputeSpecularBSDF(rough, lmodel);

    maxr = (float)abs(refl);
    if (maxr > 0.0001)
	fval += refl * specular(reflect(-vv, nn));
    return fval;
}

bsdf
gBlendTwoBSDFs(vector nn, vv;
	    vector amb1, diff1, spec1, refl1; float rough1; string lmodel1;
	    vector amb2, diff2, spec2, refl2; float rough2; string lmodel2;
	    float refl_bias; float blend_amount)
{
    bsdf	fval, aval, bval;
    bsdf	spec, aspec, bspec;
    vector	refl;
    float	maxr;
    if (lmodel1 == lmodel2)
    {
	if (lmodel1 == "constant")
	    fval = bsdf();
	else
	{
	    fval = lerp(diff1, diff2, blend_amount) * diffuse();
	    if (lmodel1 != "lambert")
	    {
		if (rough1 == rough2)
		{
		    spec = gComputeSpecularBSDF(rough1, lmodel1);
		    spec *= lerp(spec1, spec2, blend_amount);
		}
		else
		{
		    aspec = gComputeSpecularBSDF(rough1, lmodel1);
		    bspec = gComputeSpecularBSDF(rough2, lmodel2);
		    spec = (1-blend_amount)*spec1*aspec +
			      blend_amount *spec2*bspec;
		}
		fval += spec;
	    }
	}
    }
    else
    {
	// Here we're blending two separate lighting models.  So we sum up the
	// light without reflections first.
	aval = gComputeBSDF(nn, vv, amb1, diff1, spec1, {0,0,0}, rough1,
				refl_bias, lmodel1);
	bval = gComputeBSDF(nn, vv, amb2, diff2, spec2, {0,0,0}, rough2,
				refl_bias, lmodel2);
	// Until we have lerp() for BSDFs...
	//fval = lerp(aval, bval, blend_amount);
	fval = (1 - blend_amount) * aval + blend_amount * bval;
    }
    refl = lerp(refl1, refl2, blend_amount);
    maxr = (float)abs(refl);
    if (maxr > 0.0001)
	fval += refl * specular(reflect(-vv, nn));
    return fval;
}

#endif
