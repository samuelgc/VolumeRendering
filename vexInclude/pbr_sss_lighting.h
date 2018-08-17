/*
 * PROPRIETARY INFORMATION.  This software is proprietary to
 * Side Effects Software Inc., and is not to be reproduced,
 * transmitted, or disclosed in any way without written permission.
 *
 * Produced by:
 *  Side Effects Software Inc
 *  123 Front Street West
 *  Toronto, Ontario
 *  Canada   M5V 3E7
 *  416-504-9876
 *
 */
#ifndef __pbr_sss_lighting_h__
#define __pbr_sss_lighting_h__

#include "pbr.h"
#include "pbrexports.h"
#include "physicalsss.h"
#include "expsampler.h"
#include "approxsss.h"

// Returns single scatter throughput
vector
pbr_single_scatter_throughput (
    vector	inLdir; // to light
    vector	inP;
    vector	inN;
    vector	inI; // from eye
    vector	absrp;
    vector	scatr;
    float	g;
    float	eta;
    int		sid;
    float	inTime;
    string	inScope;

    // output
    vector	outP; // where light refracts into surface
    vector	outN;
    )
{
    vector dir = normalize(refract(normalize(inI), inN, eta));

    // If this is a backface, ignore
    if (dot(inN, dir) > 0.0)
	return 0.0;

    // Compute the extinction.
    vector	ext = absrp + scatr;
    float	max_dist;
    vector	pmax, nmax;

    // Trace to the next intersection.
    // TODO: This could be done by the caller instead
    if (!trace(inP, dir, inTime, "SID", sid,
		"raystyle", "nolimit",
		"scope", inScope,
		"pipeline", "displacement", 
		"P", pmax, "N", nmax))
    {
	// Simulate an infinite slab
	max_dist = 10/min(ext);
	pmax = inP + max_dist*dir;
	nmax = normalize(dir);
    }
    else
    {
	// Compute the maximum distance a ray may travel before exiting the
	// object.
	max_dist = distance(inP, pmax);
	nmax = normalize(nmax);
	if (dot(dir, nmax) < 0.0)
	    nmax = -nmax;
    }

    // Initialize the piecewise exponential sampler
    expsampler	samp;
    samp->init(ext, max_dist);

    // Initialize the radiance accumulator.
    vector radiance = 0.0;

    // Decorrelate samples from lighting.  This algorithm only uses a
    // single sample to offset the raymarch.
    float	sx, sy;
    nextsample(sid, sx, sy, "mode", "decorrelate");

    // Choose a uniform random value for depth sampling.
    float	rand = sx;
    vector	pdf;
    float	spo = samp->sample(pdf, rand);
    if (spo >= 0 && spo < max_dist)
    {
	vector psamp = inP + spo * dir;
	vector	pi, ni;

	// Trace from the sample point towards the light to
	// determine where the surface is intersected.
	if (trace(psamp, inLdir, inTime,
		  "raystyle", "nolimit",
		  "scope", inScope,
		  "bias", 0,
		  "pipeline", "displacement",
		  "P", pi,
		  "N", ni))
	{
	    ni = normalize(ni);
	}
	else
	{
	    // If we don't have a hit, we have encountered a bias
	    // error or the surface was not closed.  In these
	    // cases, we estimate the hit position/normal.
	    if (dot(inN, inLdir) >= 0.0)
	    {
		pi = inP;
		ni = normalize(inN);
	    }
	    else
	    {
		pi = pmax;
		ni = nmax;
	    }
	}

	float ldotni = abs(dot(inLdir, ni));
	float determinant = 1.0 - (eta * eta) * (1.0 - ldotni * ldotni);

	if (determinant > 0.0)
	{
	    vector	clr = 1;
	    outP = pi;
	    outN = ni;

	    // Compute the incoming light distance.
	    vector	ldir = pi - psamp;
	    float	si = length(ldir);

	    // Compute the *refracted* incoming distance.
	    si *= ldotni / sqrt(determinant);

	    // Compute the incoming fresnel factor and direction.
	    vector	ri, ti;
	    float	kri, kti;
	    fresnel(-inLdir, ni, eta, kri, kti, ri, ti);

	    // Compute the value of the phase function.  We use the
	    // Henyey-Greenstein phase function.
	    vector phase = getHenyeyGreensteinPhase(ti, -dir, g);

	    // Compute the outgoing radiance for this sample.
	    // Note that we use this sample for all color
	    // components by adjusting the sample's weight
	    // appropriately.
	    vector weight = exp(-si * ext) * phase * kti;
	    radiance += pdf * clr * weight;
	}
    }

    // Scale by the albedo.
    radiance *= scatr / ext;
    return radiance;
}

// given an array of hit distances, pick one randomly
// (weighted based on distance to the shading sample)
int _sss_select_hit(
	float pdf; // output
	float a_len[];
	float extent;
	float sx;
    )
{
    float r = sx;
    float totalWeight = 0;

    for (int i = 0; i < len(a_len); ++i)
    {
	float temp = 1.0f - abs(a_len[i] - extent) / extent;
	totalWeight += temp * temp;
    }

    r *= totalWeight;
    float w;
    int result;
    for (result = 0; result < len(a_len); ++result)
    {
	float temp = 1.0f - abs(a_len[result] - extent) / extent;
	w = temp * temp;
	r -= w;
	if (r < 0)
	    break;
    }

    pdf = w / totalWeight;
    return result;
}

//
// Generates sample point (and its throughput based on diffusion profile)
// for SSS given the mfp and sample distance parameters.
// returns 1 if successfully traced and generated SSS sample point.
// returns 0 if not.
//
int
pbr_multi_sss_throughput(
    // outputs
    vector paththroughput; 
    float colorlimit;
    vector hitN; // both input and output
    vector hitP; // both input and output
    vector hitPs;
    vector hitI;

    // inputs
    vector sssMfp;
    vector sssSampleDists;
    vector4 rnd;
    float hitTime;
    string hitScope;
    )
{
    // select a channel if the dists are not equal:
    float mfp = sssMfp.x;
    float sampleDist = sssSampleDists.x;
    if (sssSampleDists.x != sssSampleDists.y ||
	sssSampleDists.x != sssSampleDists.z)
    {
	// if channel was already chosen in previous hit
	float throughputSum = paththroughput.x + paththroughput.y + paththroughput.z;
	if (paththroughput.x == throughputSum)
	{
	    sampleDist = sssSampleDists.x;
	    mfp = sssMfp.x;
	}
	else if (paththroughput.y == throughputSum)
	{
	    sampleDist = sssSampleDists.y;
	    mfp = sssMfp.y;
	}
	else if (paththroughput.z == throughputSum)
	{
	    sampleDist = sssSampleDists.z;
	    mfp = sssMfp.z;
	}
	else
	{
	    // weigh samples based on perceptual luminance:
	    float rweight = 0.299f;
	    float gweight = 0.587f;
	    float bweight = 0.114f;

	    if (rnd.x < rweight) // red
	    {
		paththroughput.x *= 1/rweight;
		colorlimit *= gweight/rweight;
		paththroughput.y = paththroughput.z = 0;
		sampleDist = sssSampleDists.x;
	    }
	    else if (rnd.x < rweight+gweight) // green
	    {
		paththroughput.y *= 1/gweight;
		paththroughput.x = paththroughput.z = 0;
		sampleDist = sssSampleDists.y;
		mfp = sssMfp.y;
	    }
	    else // blue
	    {
		paththroughput.z *= 1/bweight;
		colorlimit *= gweight/bweight;
		paththroughput.x = paththroughput.y = 0;
		sampleDist = sssSampleDists.z;
		mfp = sssMfp.z;
	    }
	}
    }

    // fall back to simple diffuse if mfp is near zero
    if (sampleDist < PBR_MIN_LUMINANCE)
    {
	hitI = -hitN;
	return 1;
    }
    
    // generate sss sample point
    float raybias = 0;
    renderstate("renderer:raybias", raybias);

    float extent = approxsss_icdf(APPROXSSS_ICDF_N-1) * mfp + raybias;
    vector garb = {1,1,1};
    vector framen = hitN;
    vector framex, framey;
    float framenPdf, framexPdf, frameyPdf;
    makebasis(framex, framey, framen);

    // 3-way MIS
    float threewaypdf = 0;
    float sin_f = sin(rnd.y * M_PI * 2);
    float cos_f = cos(rnd.y * M_PI * 2);
    vector orig = hitP;
    vector dir;

    if (rnd.z < 0.5f)
    {
	framenPdf = 0.5f;
	framexPdf = frameyPdf = 0.25f;
    }
    else if (rnd.z < 0.75f)
    {
	vector tempv = framex;
	framex = framen;
	framen = tempv;
	framexPdf = 0.5f; 
	framenPdf = frameyPdf = 0.25f;
    }
    else
    {
	vector tempv = framey;
	framey = framen;
	framen = tempv;
	frameyPdf = 0.5f;
	framenPdf = framexPdf = 0.25f;
    }

    dir = -framen;
    orig += extent * framen + sampleDist * (framex*sin_f + framey*cos_f);
    vector prevP = hitP;

    float	a_len[];
    vector	a_n[];
    vector	a_hitPs[];
    if (trace(orig, dir, hitTime,
		"scope", hitScope,
		"maxdist", extent*2,
		"samplefilter", "all",
		"raystyle", "nolimit",
		"ray:length", a_len,
		"ray:hitPs", a_hitPs,
		"N", a_n))
    {
	float selectPdf;
	int sel = _sss_select_hit(selectPdf, a_len, extent, rnd.w);

	hitP = orig + a_len[sel] * dir;
	hitPs = a_hitPs[sel];
	hitN = normalize(a_n[sel]);

	vector toSample = normalize(hitP - prevP);
	float cavity = min(1, 1-dot(hitN, toSample));

	framenPdf *= max(abs(dot(hitN, framen)), 0.001);
	framexPdf *= abs(dot(hitN, framex));
	frameyPdf *= abs(dot(hitN, framey));
	float weight = (framenPdf * framenPdf) * cavity /
	     (framenPdf * framenPdf + framexPdf *framexPdf +  frameyPdf * frameyPdf);

	hitI = -hitN;
	float pdf = approxsss_diffusion(sampleDist, mfp) * selectPdf * framenPdf;
	paththroughput *= approxsss_diffusion(distance(prevP, hitP), mfp) * weight / pdf;
	return 1;
    }
    return 0;
}

// PBR lighting for single scatter 
void
pbr_sss_lighting (
	// exports
	vector direct;
	vector direct_comp[];
	vector direct_noshadow_comp[];

	// modified input
	float direct_samples;
	int sid; 
	pbr_sss_trace_cache ssscache;

	// inputs
	int lid;
	vector inP;
	vector inPs;
	bsdf inF;
	vector inI;
	vector inN;
	float inStep;
	float now;
	int firstbounce;

	vector tint;
	int shadow_bounces;
	int doshadow;
	float colorlimit;
	int level;
	float photonspacing;
	float samplingquality;
	string colorspace;
	int ismicropoly;
	float filter_width;
	float filter_angle;
	string currentscope;
	)
{
    int minsamples = 1;
    int maxsamples = 1;
    float threshold = 0.01;
    if (setcurrentlight(lid))
    {
	if (samplingquality > 0)
	{
	    float tmp0, tmp1, tmp2;
	    float sssqual = 1;
	    pbr_variance_settings("light:", ismicropoly,
		minsamples, maxsamples, threshold,
		tmp0, sssqual, tmp1, tmp2);
	    minsamples = max((int)(minsamples * samplingquality * sssqual), 1);
	    maxsamples = max((int)(maxsamples * samplingquality * sssqual), 1);
	}
    }

    vector ss;
    nextsample(sid, ss, "mode", "decorrelate");

    int sssseed;
    if (level == 0)
    {
	if (!ssscache.seed)
	    ssscache.seed = random_ihash(sid);
	sssseed = random_ihash(ssscache.seed);
    }
    else
    {
	sssseed = random_ihash(sid + level);
    }

    vector tmp_sum = 0;
    pbr_varianceaa var = variance_start();
    while (!var->atEnd()) 
    {
	if (var->numSamples() > 0)
	    nextsample(sid, ss, "mode", "nextpixel");

	// Sample bsdf
	int flags = PBR_BSDF_DEFAULT;
	vector eval = 0;
	float pdf = 0;
	int bsdfbounce = 0;
	vector sssSampleDists = 0;
	vector sssMfp = 0;
	vector ssAbsrp = 0;
	vector ssScatr = 0;
	float ssG = 0;
	float ssEta = 1.0;
	vector sampleDir;
	int sssType = 0;
	sample_bsdf(inF, -inI, inN, flags,
		    sampleDir, eval, pdf,
		    bsdfbounce, ss.x, ss.y, PBR_SSS_MASK,
		    // PBR SSS 
		    "import:ssstype", sssType,
		    // Multi scatter
		    "import:sssmfp", sssMfp,
		    "import:ssssampledist", sssSampleDists,
		    // Single scatter
		    "import:ssabsrp", ssAbsrp,
		    "import:ssscatr", ssScatr,
		    "import:ssg", ssG,
		    "import:sseta", ssEta
		    );

	int tmpsid = sid;
	nextsample(tmpsid, ss, "mode", "decorrelate");

	// Sample light 
	vector lpos, leval;
	float lscale;
	int light_bounces = sample_light(lid, inP, ss,
					 now, lpos, leval, lscale,
					 "N", inN,
					 "I", inI,
					 "level", level,
					 "photonspacing", photonspacing,
					 "step", inStep);
	leval *= lscale / luminance(leval);

	vector exitP, exitN;
	vector throughput = 0;
	if (sssType == PBR_SSS_TYPE_SINGLE)
	{
	    // Compute single scatter throughput
	    vector ldir = normalize(lpos - inP);
	    throughput = pbr_single_scatter_throughput(
		ldir, inP, inN, inI,
		ssAbsrp, ssScatr, ssG, ssEta,
		tmpsid, now, currentscope, exitP, exitN);
	}
	else if (sssType == PBR_SSS_TYPE_APPROX)
	{
	    int varIdx = var->numSamples();
	    if (level || len(ssscache.throughput) <= varIdx)
	    {
		// Compute multi scatter throughput:

		// NOTE that if multiple types of sss is stacked, varIdx may
		// not be contiguous and have gaps in sobol seq.
		vector4 rnd = random_sobol(2, sssseed + varIdx);
		vector tmpI;
		vector tmpthroughput = 1;
		exitN = inN;
		exitP = inP;
		vector exitPs = inPs;
		if (pbr_multi_sss_throughput(tmpthroughput, colorlimit, exitN,
		    exitP, exitPs, tmpI, sssMfp, sssSampleDists, rnd, now, 
		    currentscope))
		{
		    throughput = tmpthroughput;
		}

		if (level == 0)
		{
		    // cache
		    ssscache.exitN[varIdx]	    = exitN;
		    ssscache.exitP[varIdx]	    = exitPs;
		    ssscache.throughput[varIdx] = throughput;
		}
		exitP = exitPs; // use smoothP for shadow tests
	    }
	    else
	    {
		// fetch multiscatter throughput from cache
		exitN	    = ssscache.exitN[varIdx];
		exitP	    = ssscache.exitP[varIdx];
		throughput  = ssscache.throughput[varIdx];
	    }

	    // lambert
	    vector ldir = normalize(lpos - inP);
	    float ldotn = max(0, dot(exitN, ldir));
	    leval *= ldotn * 2.0f;
	}
	throughput *= eval;

	// Shadow
	vector noshad_comp = throughput * leval;
	vector comp = noshad_comp;
	if (luminance(noshad_comp) > PBR_MIN_LUMINANCE &&
	    doshadow && (shadow_bounces & light_bounces))
	{
	    vector shad = shadow_light(lid, exitP, lpos-inP, now,
		"nofakecaustics", 0,
		"noantialiasing", 1,
		"N", exitN,
		"SID", tmpsid,
		"filter_width", filter_width,
		"filter_angle", filter_angle);
	    comp *= shad;
	}

	FOR_ALL_EXPORTS
	{
	    int bmask = bouncemask(cname);
	    if ( (!firstbounce && bsdfbounce & bmask) ||
		 (firstbounce & bmask) )
	    {
		direct_comp[cidx] += tint * comp;
		direct_noshadow_comp[cidx] += tint * noshad_comp;
	    }
	}

	tmp_sum += tint * comp;
	var->advance(PBR_ALL_MASK, tmp_sum, minsamples, maxsamples, colorspace, threshold, 1, 1, 1, 1);
    }

    if (var->numSamples())
    {
	float isamples = 1.0 / (float)var->numSamples();
	direct = tmp_sum * isamples;
	direct_samples += var->numSamples();

	for (int i = 0; i < len(direct_comp); i++)
	{
	    direct_comp[i] *= isamples;
	    direct_noshadow_comp[i] *= isamples;
	}
    }
}

#endif
