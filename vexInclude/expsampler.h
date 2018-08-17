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
 * NAME:	expsampler.h ( VEX )
 *
 * COMMENTS:	A high quality sampling algorithm for a 3-function
 *		piecewise exponential (usually RGB attenuation)
 */

#ifndef __expsampler__
#define __expsampler__

// Sort vector components in descending order.  So (1,2,3) would become
// (3,2,1).
vector
sortDescending(const vector val)
{
    float	mid;

    if (val[0] > val[1])
	mid = (val[2] > val[0] ? val[0] : max(val[1], val[2]));
    else
	mid = (val[2] > val[1] ? val[1] : max(val[0], val[2]));

    return set(max(val), mid, min(val));
}

// Solve for the crossing points.  Given 2 exponentials:
//	a*exp(-a*x)
//	b*exp(-b*x)
// The intersection is:
//	x = log(a/b)/(a-b)
float
intersectExp(const float a; const float b; const float def)
{
    return a > b ? log(a/b)/(a-b) : def;
}

// Return the sampled component for a 3-valued cdf
int
sampleExpComp(const float val; const float k1; const float k2)
{
    return select(val < k1, 0, select(val < k2, 1, 2));
}

struct expsampler
{

    // Create a new exponential sampler for the 3 exponents defined by the ext
    // vector:
    //	exp(-ext*x)
    void init(const vector inext; const float maxdist)
    {
	ext = max(inext, 1e-6F);

	// Sort vector components in descending order.  We'll construct a
	// piecewise continuous pdf from the 3 exponential functions by finding
	// the crossing points (x1 and x2) between the different pdfs.
	ext_sort = sortDescending(ext);

	x1 = intersectExp(ext_sort[0], ext_sort[1], 0);
	x2 = intersectExp(ext_sort[1], ext_sort[2], x1);

	// Solve for integrals in the 3 ranges.  The resulting integral from 0
	// to infinity (k3) will be used as the normalization factor for the
	// pdf.

	// Integral from 0..x1, 0..x2, 0..infinity
	k1 = 1 - exp(-ext_sort[0]*x1);
	k2 = exp(-ext_sort[1]*x1) -
		  exp(-ext_sort[1]*x2) + k1;
	k3 = exp(-ext_sort[2]*x2) + k2;

	// Integration constant at 0, x1, x2
	o[0] = 1;
	o[1] = k1 + exp(-ext_sort[1]*x1);
	o[2] = k3;

	// Compute the maximum uniform random value that will result in a
	// sample distance that lies inside the model.
	int icomp = sampleExpComp(maxdist, x1, x2);
	max_rand = (o[icomp] -
		exp(-ext_sort[icomp]*maxdist)) / k3;
    }

    // Distribute a sample point along the refracted outgoing ray according
    // to the piecewise exponential function:
    //	pdf: exp(-ext_sort[icomp] * dist) * ext_sort[icomp] / k3
    //	cdf: (1 - exp(-ext_sort[icomp] * dist)) / k3.
    // For:
    //	icomp := 0	:       x < x1
    //	icomp := 1	: x1 <= x < x2
    //	icomp := 2	: x2 <= x
    float sample(export vector pdf; const float rand)
    {
	float sval = rand * k3 * max_rand;

	int icomp = sampleExpComp(sval, k1, k2);
	float spo = -log(o[icomp] - sval) / ext_sort[icomp];

	// Calculate the pdf for each component
	pdf = exp(-ext*spo) * ext;

	// Adjust weight due to importance sampling
	pdf *= k3 * max_rand;
	pdf /= exp(-ext_sort[icomp]*spo) * ext_sort[icomp];

	return spo;
    }

    vector	ext;		// Exponent
    vector	ext_sort;	// Sorted exponent
    float	max_rand;	// Maximum random value
    float	x1, x2;		// Crossing distances
    float	k1, k2, k3;	// Integrals
    vector	o;		// Offsets
};

#endif
