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
 * NAME:	pbr_direct_lighting.h ( VOP Library Functions, VEX )
 *
 * COMMENTS:	This file provides a function to accumulate lighting over
 *		all lights.
 */

#ifndef __pbr_direct_lighting_h__
#define __pbr_direct_lighting_h__

#include "pbr.h"
#include "pbr_sss_lighting.h"
#include "pbrexports.h"

///
/// pbr_direct_lighting() will compute direct lighting and sum the results into
/// the output arrays:
///	direct_light
///	direct_{exports}[]
///
/// Direct lighting is computed by sampling lights and the bsdf and tracing
/// up to 2 shadow rays - one for a light-centric sample, and another for a
/// surface-centric sample.  Multiple importance sampling is used to weight
/// the samples.
///
/// This method iterates over the lights, calling the illuminate() shader (as
/// defined by the vm_illumshader property).
///
/// Input parameters used for light selection:
///  - material inMat, vector inP @n
///	The light mask, light categories of the @c inMat parameter are used for
///	rejection testing of the lights (combined with the active radius as
///	determined by @c inP).
///
/// Input parameters passed to light->illuminate() shader.
///  - int inSID @n
///	Sample ID associated with this lighting computation
///  - vector inP @n
///	Point on the surface being lit
///  - bsdf inF @n
///	BSDF for the surface being lit
///  - vector inI, inN @n
///	Incident and normal vectors for the surface being lit
///  - float inStep @n
///	TODO
///  - float now @n
///	TODO
///  - vector tint
///
void
pbr_direct_lighting(
	// Exports
	vector direct_light;	// Shadowed light
	vector direct_comp[];	// Shadowed light (per-component)
	vector direct_noshadow_comp[];	// Unshadowed light (per-component)
	float direct_samples;
	int indirect_bounces;	// Indirect bounces to mask out
	pbr_sss_trace_cache ssscache;

	// Globals
	material inMat;
	int inSID;
	vector inP;
	vector inPs;
	bsdf inF;
	vector inI;
	vector inN;
	float inStep;
	float now;

	// Inputs
	vector tint;
	int lights[];
	int light_map[];	// Map from lid to export array index
	int inbounces;
	int shadow_bounces;
	int firstbounce;
	pbr_lighting_parms lparms;
	float colorlimit;
	int level;
	float photonspacing;
	float samplingquality;
	string colorspace;
	int ismicropoly;
	float rayfilteramount;
	float filter_width;
	float filter_angle;
	float maxrough)
{
    // SSS should not be handled by illum shader
    int		bounces = inbounces & ~PBR_SSS_MASK;
    int		sid = inSID;

    indirect_bounces = 0;

    int nlights = len(lights);

    int
    shouldProcessLight(material inMat; vector inP; int lid; int level)
    {
	if (level == 0)
	    return haslight(inMat, inP, lid, "direct", 1);
	return haslight(inMat, inP, lid, "direct", 0);
    }

    int hasSSS = luminance(albedo(inF, -inI, inbounces & PBR_SSS_MASK, "direct", 1)) > 
	PBR_MIN_LUMINANCE;
    string currentHitScope = "";
    if (hasSSS)
	renderstate(inMat, "object:name", currentHitScope);
    
    for (int i = 0; i < nlights; i++)
    {
	int	lid = lights[i];

	if (!shouldProcessLight(inMat, inP, lid, level))
	    continue;

	// Evaluate into write-only temporaries
	vector	tmp_light;
	vector	tmp_noshadow_comp[];
	float	tmp_samples;
	vector	tmp_comp[];

	light	lt = getlight(lid);

	// Single scatter case
	if (hasSSS)
	{
	    vector direct = 0;
	    float tmpcolorlimit = colorlimit;
	    pbr_sss_lighting (direct, tmp_comp, tmp_noshadow_comp, 
		direct_samples, sid, ssscache, lid, inP, inPs, inF, inI, inN,
		inStep, now, firstbounce, tint, shadow_bounces, lparms.doshadow,
		tmpcolorlimit, level, photonspacing, samplingquality,
		colorspace, ismicropoly, filter_width, filter_angle,
		currentHitScope);

	    direct_light += direct;
	    int lmapidx = light_map[lid];
	    FOR_ALL_EXPORTS
	    {
		LCOMP(direct_comp, lmapidx, cidx) += tmp_comp[cidx];
		LCOMP(direct_noshadow_comp, lmapidx, cidx) += tmp_noshadow_comp[cidx];
	    }
	}

	if (bounces)
	{
	    lt->illuminate(
		    "lid", lid,

		    "direct_light", tmp_light,
		    "direct_noshadow_comp", tmp_noshadow_comp,
		    "direct_samples", tmp_samples,
		    "direct_comp", tmp_comp,

		    "bounces", bounces,
		    "indirect_bounces", indirect_bounces,
		    "sid", sid,

		    "inP", inP,
		    "inPs", inPs,
		    "inF", inF,
		    "inI", inI,
		    "inN", inN,
		    "inStep", inStep,
		    "now", now,
		    "tint", tint,
		    "shadow_bounces", shadow_bounces,
		    "firstbounce", firstbounce,
		    "doshadow", lparms.doshadow,
		    "colorlimit", colorlimit,
		    "level", level,
		    "photonspacing", photonspacing,
		    "samplingquality", samplingquality,
		    "colorspace", colorspace,
		    "ismicropoly", ismicropoly,
		    "rayfilteramount", rayfilteramount,
		    "filter_width", filter_width,
		    "filter_angle", filter_angle,
		    "maxrough", maxrough
		    );

	    int lmapidx = light_map[lid];

	    direct_light += tmp_light;
	    direct_samples += tmp_samples;
	    FOR_ALL_EXPORTS
	    {
		LCOMP(direct_comp, lmapidx, cidx) += tmp_comp[cidx];
		LCOMP(direct_noshadow_comp, lmapidx, cidx) += tmp_noshadow_comp[cidx];
	    }
	}

	if (!bounces && !hasSSS)
	    break;
    }
}

#endif
// :vimvfl:
