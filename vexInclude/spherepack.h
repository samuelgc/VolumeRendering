/*
 * PROPRIETARY INFORMATION.  This software is proprietary to
 * Side Effects Software Inc., and is not to be reproduced,
 * transmitted, or disclosed in any way without written permission.
 *
 * Produced by:
 *  Side Effects Software Inc
 *  477 Richmond Street West
 *  Toronto, Ontario
 *  Canada   M5V 3E7
 *  416-504-9876
 *
 * NAME:    spherepack.h ( VOP Library Functions, VEX )
 *
 * COMMENTS:    Support functions for PointsFromVolume SOP.
 */

#ifndef __spherepack_h__
#define __spherepack_h__

#include <math.h>

// Returns 1 if we can skip point generation for this bounding box
// given the jitter, SDF, and iso values.
int skipbounds(vector bbmin; vector bbmax; vector jitter;
               int vol; string volname; float miniso; float maxiso)
{
    vector size = bbmax - bbmin;
    vector center = bbmin + size / 2;
    float r = length(size) / 2;
    float maxdist = r + max(jitter);
    float dist = volumesample(vol, volname, center);
    return (dist > maxdist + maxiso || dist < -maxdist + miniso);
}

// Point jitter.
// TODO - make more robust for huge point counts.
vector pointjitter(vector p; vector jitter)
{
    return fit01(rand(p), -jitter, jitter);
}

int jitteredlattice(int geo; vector psep; vector bbmin; vector bbmax; 
                 vector jitter; float jitterseed; vector offset; float reject;
                 int vol; string volname; float miniso; float maxiso;
		 int dodither; vector dithernormal; float ditherangle; float dithersep)
{
    int npts = 0;

    // Add some quantization to account for possible rounding errors on the
    // input bboxes. For example if all the bboxes come from voxels in a
    // Volume VOP, the boundaries might fall on the wrong side of roundoff
    // error for the multiple adjacent boxes, particularly when voxel size is
    // a multiple of particle separation, as is often the case.
    vector quantize = 0.001 * psep;
    vector begin = ceil((bbmin - offset + quantize) / psep);
    vector end = floor((bbmax - offset + quantize) / psep);
    int dojitter = length(jitter) > 0;
    int doreject = reject > 0;
    vector vecseed = set(jitterseed, 
                        (jitterseed + 123) * M_SQRT_2_3, 
                        (jitterseed + 457) * M_2SQRT6_3);
    for(int i=(int)begin.x; i <= end.x; i++)
    for(int j=(int)begin.y; j <= end.y; j++)
    for(int k=(int)begin.z; k <= end.z; k++)
    {
        vector idx = set(i, j, k);
        // Probabilistically reject samples if requested.
        if (doreject && rand(idx + vecseed * M_SQRT3) < reject)
            continue;
        vector samplepos = psep * idx + offset;
        if (dojitter)
            samplepos += pointjitter(idx + vecseed, jitter);
        // Create point if inside SDF depth limits.
        float sampledist = volumesample(vol, volname, samplepos);
	if (dodither)
	{
	    int		keep = 0;
	    if (sampledist <= maxiso - dithersep && sampledist >= miniso + dithersep)
	    {
		// Guaranteed in...
		keep = 1;
	    }
	    else if (sampledist <= maxiso + dithersep && sampledist >= miniso - dithersep)
	    {
		// Potential dither

		int		shoulddither = 1;
		if (ditherangle > -1)
		{
		    // Determine if our local normal is inside the
		    // dither threshold.
		    vector vnormal = volumegradient(vol, volname, samplepos);
		    vnormal = normalize(vnormal);
		    shoulddither = 0;
		    if (dot(vnormal, dithernormal) >= ditherangle)
		    {
			shoulddither = 1;
		    }
		}

		if (shoulddither)
		{
		    // Determine if we are kept by random rejection
		    float error = max(sampledist - (maxiso-dithersep),
				      miniso + dithersep - sampledist);
		    error /= (dithersep * 2);
		    error = 1 - error;

		    if (rand(vecseed + samplepos) < error)
			keep = 1;
		}
		else
		{
		    // Use hard edge condition.
		    if (sampledist <= maxiso && sampledist >= miniso)
		    {
			keep = 1;
		    }
		}
	    }
	    else
	    {
		// Guaranteed out.
	    }

	    if (keep)
	    {
		addpoint(geo, samplepos);
		npts++;
	    }
	}
	else
	{
	    if (sampledist <= maxiso && sampledist >= miniso)
	    {
		addpoint(geo, samplepos);
		npts++;
	    }
	}
    }
    return npts;
}

int simplecubicpts(int geo; vector psep; vector bbmin; vector bbmax; 
                 vector jitter; float jitterseed; vector offset; float reject;
                 int vol; string volname; float miniso; float maxiso;
		 int dodither; vector dithernormal; float ditherangle; float dithersep)
{
    if (skipbounds(bbmin, bbmax, jitter, vol, volname, miniso, maxiso))
        return 0;
    return jitteredlattice(geo, psep, bbmin, bbmax,
                               jitter, jitterseed, offset, reject,
                               vol, volname, miniso, maxiso,
			       dodither, dithernormal, ditherangle, dithersep);
}

int tetrahedralpts(int geo; vector psep; vector bbmin; vector bbmax; 
                 vector jitter; float jitterseed; vector offset; float reject;
                 int vol; string volname; float miniso; float maxiso;
		 int dodither; vector dithernormal; float ditherangle; float dithersep)
{
    if (skipbounds(bbmin, bbmax, jitter, vol, volname, miniso, maxiso))
        return 0;
    int npts = 0;
    vector localoff = offset;
    vector localpsep = psep * { 1, M_2SQRT6_3, M_SQRT3 };

    // layer 1
    npts += jitteredlattice(geo, localpsep, bbmin, bbmax,
                    jitter, jitterseed, localoff, reject,
                    vol, volname, miniso, maxiso,
		    dodither, dithernormal, ditherangle, dithersep);
    // layer 2
    localoff += psep / 2 * { 1, 0, M_SQRT3 };
    npts += jitteredlattice(geo, localpsep, bbmin, bbmax, 
                    jitter, (jitterseed + 1) * M_PI, localoff, reject,
                    vol, volname, miniso, maxiso,
		    dodither, dithernormal, ditherangle, dithersep);
    // layer 4
    localoff += psep * { 0, M_SQRT_2_3,  M_1_SQRT3};
    npts += jitteredlattice(geo, localpsep, bbmin, bbmax,
                    jitter, (jitterseed + 2) * M_E, localoff, reject,
                    vol, volname, miniso, maxiso,
		    dodither, dithernormal, ditherangle, dithersep);
    // layer 3
    localoff = offset + psep * { 0, M_SQRT_2_3,  M_1_SQRT3};
    npts += jitteredlattice(geo, localpsep, bbmin, bbmax,
                    jitter, (jitterseed + 3) * M_SQRT2, localoff, reject,
                    vol, volname, miniso, maxiso,
		    dodither, dithernormal, ditherangle, dithersep);
    return npts;
}


int spherepackbbox(int geo; float psep; vector bbmin; vector bbmax; int pattern;
                vector jitter; float jitterseed; vector offset;
                int vol; string volname; float miniso; float maxiso;
                float oversample; float oversampledist;
		int dodither; vector dithernormal; float ditherangle)
{
    int npts = 0;
    if (pattern == 0)
    {
        npts += simplecubicpts(geo, psep, bbmin, bbmax,
                                   jitter, jitterseed, offset, 0,
                                   vol, volname, miniso, maxiso,
				   dodither, dithernormal, ditherangle, psep);
    }
    else if (pattern == 1)
    {
        npts += tetrahedralpts(geo, psep, bbmin, bbmax,
                                   jitter, jitterseed, offset, 0,
                                   vol, volname, miniso, maxiso,
				   dodither, dithernormal, ditherangle, psep);
    }
    // Oversample if necessary.
    oversample = max(oversample, 0);
    int noversample = (int)ceil(oversample);
    // How many points to reject as we're looping to oversample.
    float reject = 1 - (oversample / noversample);
    // Clamp oversample miniso to regular miniso.
    float oversampleminiso = max(maxiso - oversampledist, miniso);
    float seed = 347.137;
    for (int i=0; i < noversample; i++)
    {
        if (pattern == 0)
        {
            npts += simplecubicpts(geo, psep, bbmin, bbmax,
                                   jitter, jitterseed + seed, offset, reject,
                                   vol, volname, oversampleminiso, maxiso,
				   dodither, dithernormal, ditherangle, psep);
        }
        else if (pattern == 1)
        {
            npts += tetrahedralpts(geo, psep, bbmin, bbmax,
                                   jitter, jitterseed + seed, offset, reject,
                                   vol, volname, oversampleminiso, maxiso,
				   dodither, dithernormal, ditherangle, psep);
        }
        seed *= M_E;
    }
    return npts;
}

#endif

