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
 * NAME:	pbr.h ( VOP Library Functions, VEX )
 *
 * COMMENTS:	Support functions for PBR.  Most of these library functions
 *		are shading context only.
 */

#ifndef __pbr_h__
#define __pbr_h__

#include "math.h"

//-----------------------------------------------------------------------------
// Bounce masks

// Ambient mask is only used to identify ambient lights when performing

// direct lighting
#define PBR_AMBIENT_MASK	bouncemask("ambient")

#define PBR_DIFFUSE_MASK	bouncemask("alldiffuse")
#define PBR_REFRACT_MASK	bouncemask("allrefract")
#define PBR_VOLUME_MASK		bouncemask("allvolume")
#define PBR_SSS_MASK            bouncemask("allsss")

#define PBR_NO_MASK		0
#define PBR_ALL_MASK		bouncemask("all")

#define PBR_GLOSSY_MASK		\
    (PBR_ALL_MASK & ~(PBR_DIFFUSE_MASK|PBR_VOLUME_MASK|PBR_SSS_MASK))

// Reflection bounce types are defined as all other bounce types, to ensure
// we have full coverage of all user-defined components.
#define PBR_REFLECT_MASK		\
    (PBR_GLOSSY_MASK & ~PBR_REFRACT_MASK)

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// BSDF split types

#define PBR_SPLIT_FULL          0
#define PBR_SPLIT_RANDOM        1
#define PBR_SPLIT_ALBEDO        2
#define PBR_SPLIT_COMPONENT     3
#define PBR_SPLIT_DEFAULT       PBR_SPLIT_ALBEDO

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// BSDF flags

#define PBR_BSDF_I_EVENT_ENTER      0x0001
#define PBR_BSDF_I_EVENT_INTERIOR   0x0002
#define PBR_BSDF_I_EVENT_EXIT       0x0004
#define PBR_BSDF_O_EVENT_ENTER      0x0008
#define PBR_BSDF_O_EVENT_INTERIOR   0x0010
#define PBR_BSDF_O_EVENT_EXIT       0x0020
#define PBR_BSDF_DELTA              0x0040
#define PBR_BSDF_DIRECT             0x0080
#define PBR_BSDF_REVERSE            0x0100
#define PBR_BSDF_MONO               0x0200

// Note: If you intend to extend the list of BSDF flags for your use
// please ensure the following flags are not used.
#define PBR_BSDF_RESERVED_5         0x0400
#define PBR_BSDF_RESERVED_4         0x0800
#define PBR_BSDF_RESERVED_3         0x1000
#define PBR_BSDF_RESERVED_2         0x2000
#define PBR_BSDF_RESERVED_1         0x4000
#define PBR_BSDF_RESERVED_0         0x8000

#define PBR_BSDF_DEFAULT            0

#define PBR_BSDF_I_EVENT_ANY \
     (PBR_BSDF_I_EVENT_ENTER \
    | PBR_BSDF_I_EVENT_INTERIOR \
    | PBR_BSDF_I_EVENT_EXIT) \

#define PBR_BSDF_O_EVENT_ANY \
     (PBR_BSDF_O_EVENT_ENTER \
    | PBR_BSDF_O_EVENT_INTERIOR \
    | PBR_BSDF_O_EVENT_EXIT) \

// SSS type
#define PBR_SSS_TYPE_NONE	    0
#define PBR_SSS_TYPE_APPROX	    1
#define PBR_SSS_TYPE_SINGLE	    2

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Primitive shading style types

#define PBR_PRIMSS_NONE         0 // Can't shade
#define PBR_PRIMSS_GRID         1 // Rectangular grid
#define PBR_PRIMSS_POINTS       2 // Points
#define PBR_PRIMSS_PSOUP        3 // Polygon soups (metaballs)
#define PBR_PRIMSS_VOLUME       4 // Volume

//-----------------------------------------------------------------------------

// Min/max luminance values
#define PBR_MIN_LUMINANCE	1e-6
#define PBR_MAX_LUMINANCE	1e30

// Get index according to masks set in 'PATHSTATE'.
#define PBR_PATHSTATE_IDX(PATHSTATE) \
    ( PATHSTATE & (PBR_DIFFUSE_MASK | PBR_VOLUME_MASK) ) ? 0 \
  : ( PATHSTATE & PBR_SSS_MASK )                         ? 1 \
  : ( PATHSTATE & PBR_REFRACT_MASK )                     ? 2 \
  : ( PATHSTATE & PBR_ALL_MASK )                         ? 3 \
  : 0

//
// Return the pixel coordinates for a given point in arbitrary space.
//

vector
pbr_get_pixel(vector p;)
{
    vector res;
    renderstate("renderer:resolution", res);

    vector pr = ptransform("space:current", "space:ndc", p);
    pr.x *= res.x;
    pr.y *= res.y;
    pr.x = floor(pr.x);
    pr.y = floor(pr.y);
    pr.z = 0;

    return pr;
}

//
// Return a sample uniformly distributed on a disc of arbitrary
// radius, with falloff from the center.
//

vector
pbr_sample_disc(float u;
                float v;
                float radius;
                float radiuspower)
{
    float r = sqrt(pow(u, radiuspower)) * radius;
    float theta = M_TWO_PI * v;
    vector x;
    x.x = r * cos(theta);
    x.y = r * sin(theta);
    x.z = 0;

    return x;
}

void
pbr_increment_level(
	int pathstate;
	int reflectlevel;
	int refractlevel;
	int diffuselevel;
        int ssslevel;
	int volumelevel;
	string raystyle)
{
    pathstate &= PBR_ALL_MASK;
    if (pathstate & PBR_DIFFUSE_MASK)
    {
	diffuselevel++;
	raystyle = "diffuse";
    }
    else if (pathstate & PBR_SSS_MASK)
    {
        ssslevel++;
        raystyle = "sss";
    }
    else if (pathstate & PBR_REFLECT_MASK)
    {
	reflectlevel++;
	raystyle = "reflect";
    }
    else if (pathstate & PBR_REFRACT_MASK)
    {
	refractlevel++;
	raystyle = "refract";
    }
    else if (pathstate & PBR_VOLUME_MASK)
    {
	volumelevel++;
	raystyle = "diffuse";
    }
}

int
pbr_findlight(int lights[]; int lid; int end)
{
    for (int i = 0; i < end; i++)
	if (lights[i] == lid)
	    return i;
    return -1;
}

float
pbr_filterangle(float pdf; float amount)
{
    return amount / sqrt(pdf*16 + 4);
}

void
pbr_variance_settings(string prefix; int ismicropoly;
              int minsamples; int maxsamples;
              float threshold;
              float diffusequality;
              float sssquality;
              float refractionquality;
              float reflectionquality;)
{
    minsamples = maxsamples = 1;
    renderstate(prefix + "minraysamples", minsamples);
    renderstate(prefix + "maxraysamples", maxsamples);

    threshold = 0;
    if (maxsamples > minsamples)
    {
        renderstate(prefix + "variance", threshold);
        renderstate(prefix + "diffusequality", diffusequality);
        renderstate(prefix + "sssquality", sssquality);
        renderstate(prefix + "refractionquality", refractionquality);
        renderstate(prefix + "reflectionquality", reflectionquality);
    }

    if (!ismicropoly)
    {
        vector	vsamples = 0;
        renderstate("renderer:samples", vsamples);
        float	psamples = vsamples.x * vsamples.y;
        threshold *= sqrt(psamples);
    }
}

void
pbr_indirect_variance_settings(string prefix; int ismicropoly;
		      int minsamples; int maxsamples;
                      float threshold;
                      float diffusequality;
                      float sssquality;
                      float refractionquality;
                      float reflectionquality;
                      float volumequality;)
{
    minsamples = maxsamples = 1;
    int decoupleindirect = 0;
    renderstate(prefix + "decoupleindirect", decoupleindirect);
    if (decoupleindirect)
    {
        renderstate(prefix + "minindirectraysamples", minsamples);
        renderstate(prefix + "maxindirectraysamples", maxsamples);
    }
    else
    {
        renderstate(prefix + "minraysamples", minsamples);
        renderstate(prefix + "maxraysamples", maxsamples);
    }

    threshold = 0;
    float globalquality = 1.0;
    diffusequality = 1.0;
    sssquality = 1.0;
    refractionquality = 1.0;
    reflectionquality = 1.0;
    if (maxsamples > minsamples)
    {
        if (decoupleindirect)
            renderstate(prefix + "indirectvariance", threshold);
        else
            renderstate(prefix + "variance", threshold);

        renderstate(prefix + "globalquality", globalquality);
        renderstate(prefix + "diffusequality", diffusequality);
        renderstate(prefix + "sssquality", sssquality);
        renderstate(prefix + "refractionquality", refractionquality);
        renderstate(prefix + "reflectionquality", reflectionquality);
        renderstate(prefix + "volumequality", volumequality);

        diffusequality *= globalquality;
        sssquality *= globalquality;
        refractionquality *= globalquality;
        reflectionquality *= globalquality;
        volumequality *= globalquality;
    }

    if (!ismicropoly)
    {
	vector	vsamples = 0;
	renderstate("renderer:samples", vsamples);
	float	psamples = vsamples.x * vsamples.y;
	threshold *= sqrt(psamples);
    }
}

//
// Grid discretizing path tracing data (float variant)
//

struct pbr_grid_f
{
    void resize(int resolution_;)
    {
        resolution = max(1, resolution_);
        dx = 1.0 / resolution;
        resize(buf, resolution*resolution);
    }

    void assign(float value;)
    {
        for (int i = 0; i < resolution*resolution; i++)
        {
            buf[i] = value;
        }
    }

    void addSample(float sx;
                   float sy;
                   float value;)
    {
        int ix = floor(sx * resolution);
        int iy = floor(sy * resolution);
        buf[iy*resolution + ix] += value;
    }

    void scale(float x;)
    {
        for (int i = 0; i < resolution*resolution; i++)
        {
            buf[i] *= x;
        }
    }

    float sum()
    {
        float x = 0;

        for (int i = 0; i < resolution*resolution; i++)
        {
            x += buf[i];
        }

        return x;
    }

    float average()
    {
        float x = this->sum();

        return x / (resolution*resolution);
    }

    void range(float vmin, vmax)
    {
        vmin = 1e6;
        vmax = -1e6;

        for (int i = 0; i < resolution*resolution; i++)
        {
            vmin = min(vmin, buf[i]);
            vmax = max(vmax, buf[i]);
        }
    }

    float variance()
    {
        float vmin = 1e6;
        float vmax = -1e6;
        this->range(vmin, vmax);

        return vmax - vmin;
    }

    void write(string filename;)
    {
        int i = 0;

        for (int y = 0; y < resolution; y++)
        {
            for (int x = 0; x < resolution; x++)
            {
                vector Pr;
                Pr.x = x;
                Pr.y = y;
                Pr.z = 0;

                vector Ce;
                Ce.x = buf[i];
                Ce.y = buf[i];
                Ce.z = buf[i];

                pcwrite(filename, "P", Pr, "Ce", Ce);

                i++;
            }
        }
    }

    int resolution;
    float dx;
    float buf[];
}

//
// Object to simplify variance antialiasing
//

struct pbr_varianceaa
{
    void reset()
    {
        var = { 0, 0, 0, 0 };
        prevlum = { 0, 0, 0, 0 };
        done = 0;
        nsamples = { 0, 0, 0, 0 };
    }

    int numSamples()
    {
        return (nsamples[0] + nsamples[1] + nsamples[2] + nsamples[3]);
    }

    int atEnd()
    {
        return done;
    }

    void advance(int pathstate; vector value;
         int minsamples; int maxsamples;
         string colorspace;
         float threshold;
         float diffusequality;
         float sssquality;
         float refractionquality;
         float reflectionquality;)
    {
        int psi = PBR_PATHSTATE_IDX(pathstate);

        int wmaxsamples[] = {0, 0, 0, 0};
        wmaxsamples[0] = floor(maxsamples * diffusequality);
        wmaxsamples[1] = floor(maxsamples * sssquality);
        wmaxsamples[2] = floor(maxsamples * refractionquality);
        wmaxsamples[3] = floor(maxsamples * reflectionquality);

        if ((wmaxsamples[psi] > minsamples) && (nsamples[psi]+1) < wmaxsamples[psi])
        {
            float lum = luminance(value) / (nsamples[psi]+1);
            if (colorspace == "gamma")
            lum = sqrt(lum);

            int		samplesize;
            float	mean;
            float	newvar = variance(lum - prevlum[psi], mean, samplesize);

            var[psi] = (var[psi]*nsamples[psi] + newvar) / (nsamples[psi]+1);
            prevlum[psi] = lum;

            done = (nsamples[psi]+1) >= minsamples
                    && var[0] < (threshold*threshold * 1.0 / diffusequality)
                    && var[1] < (threshold*threshold * 1.0 / sssquality)
                    && var[2] < (threshold*threshold * 1.0 / refractionquality)
                    && var[3] < (threshold*threshold * 1.0 / reflectionquality);
        }
        else
        {
            done = (nsamples[psi]+1) >= minsamples;
        }

        nsamples[psi]++;
    }

    float   var[] = { 0, 0, 0, 0 };
    float   prevlum[] = { 0, 0, 0, 0 };
    int     done;
    int     nsamples[] = { 0, 0, 0, 0 };
}

pbr_varianceaa variance_start()
{
    return pbr_varianceaa();
}

// Parameters to control how lighting is computed, grouped together in a
// struct to make it easier to pass around
struct pbr_lighting_parms {
    int		doshadow;
}

//
// Calculate allowable bounces at given pathtracing levels for a given type
// of path
//

void
pbr_bounce_mask(int direct_bounces; int indirect_bounces;
		int shadow_bounces; int background_bounces;
		int direct_mask; int indirect_mask;
		string pathtype; int pathstate; int hitBounces;
		int reflectlevel; int reflectlimit;
		int refractlevel; int refractlimit;
		int diffuselevel; int diffuselimit;
                int ssslevel; int ssslimit;
		int volumelevel; int volumelimit;
		string raylimiteval; int raylimitmask;)
{
    // Find out what paths we need to trace
    if (pathtype == "diffuse" ||
	pathtype == "specular")
    {
        if (pathstate & (PBR_DIFFUSE_MASK | PBR_VOLUME_MASK | PBR_SSS_MASK))
	{
	    // Once we have one diffuse bounce, need to
	    // continue sampling only diffuse bounces.  If
	    // we've hit the diffuse bounce limit, then make
	    // sure that we only compute direct lighting.
            indirect_bounces = PBR_DIFFUSE_MASK | PBR_VOLUME_MASK | PBR_SSS_MASK;
            direct_bounces = PBR_DIFFUSE_MASK | PBR_VOLUME_MASK | PBR_SSS_MASK;
	}
	else
	{
	    // If we're using diffuse bounces, the first will
	    // sample any BSDF component.  Subsequent bounces
	    // will only sample diffuse or specular, depending
	    // on the first bounce type unless we have
	    // caustics.
	    indirect_bounces = PBR_ALL_MASK;
	    direct_bounces = PBR_ALL_MASK;
	}
    }
    else
    {
	indirect_bounces = PBR_ALL_MASK;
	direct_bounces = PBR_ALL_MASK;
    }

    // Mask out disallowed bounces
    direct_bounces &= direct_mask;
    indirect_bounces &= indirect_mask;

    // Mask out specific bounces if we've exceeded the bounce limit
    raylimitmask = 0;
    if (reflectlevel >= reflectlimit)
        raylimitmask |= PBR_REFLECT_MASK;
    if (refractlevel >= refractlimit)
        raylimitmask |= PBR_REFRACT_MASK;
    if (diffuselevel >= diffuselimit)
        raylimitmask |= PBR_DIFFUSE_MASK;
    if (ssslevel >= ssslimit)
    {
        raylimitmask |= PBR_SSS_MASK;
        //direct_bounces &= ~PBR_SSS_MASK;
    }
    if (volumelevel >= volumelimit)
        raylimitmask |= PBR_VOLUME_MASK;

    background_bounces = 0;
    shadow_bounces = direct_bounces;
    if (raylimiteval == "direct")
    {
        indirect_bounces &= ~(raylimitmask & ~PBR_GLOSSY_MASK);
        raylimitmask &= PBR_GLOSSY_MASK;
        shadow_bounces &= ~raylimitmask;
        background_bounces = raylimitmask;
    }
    else
    {
        indirect_bounces &= ~raylimitmask;
    }

    // Mask out bounces that are not part of the hit bsdf
    direct_bounces &= hitBounces;
    indirect_bounces &= hitBounces;
}

float
pbr_shadowmattecomp(float lit; float shad)
{
    return shad > 1e-6 ? max(1.0 - lit/shad, 0.0) : 0.0;
}

vector
pbr_shadowmatte(vector lit; vector shad)
{
    return (vector)pbr_shadowmattecomp(max(lit), max(shad));
}

// Clamp k and eval to the given color limit.  Here, k represents the value
// of a pdf.
float
pbr_clampCf(vector clr; float colorlimit)
{
    float	kscale = 1;
    float	lum = luminance(clr);
    if (colorlimit >= 0 && lum > colorlimit)
	kscale = colorlimit/lum;
    return kscale;
}

struct pbr_sss_trace_cache
{
    int seed = 0;
    vector exitN[];
    vector exitP[];
    vector throughput[];
};

#endif
// :vimvfl:
