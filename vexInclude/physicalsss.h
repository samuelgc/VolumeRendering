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
 * NAME:	physicalsss.h ( VEX )
 *
 * COMMENTS:	Functions for physically based sub-surface scattering.
 *
 * Many functions share a common set of parameters.  Rather than documenting
 * each parameter repeatedly in each function comment block, we describe those
 * parameters here:
 *
 * absrp:
 * 	Absorption coefficient that specifies the probability per unit distance
 * 	that an absorption event will occur.  Since the probability is per unit
 * 	distance, this parameter has a range of [0, +inf].  Higher values
 * 	correspond to more absorption.  This parameter is a vector and each
 * 	component corresponds to a different wavelength (R, G, B).
 *
 * albedo:
 * 	The ratio of scattering to absorption plus scattering:
 *
 * 		albedo = scatr / (scatr + absrp)
 *
 * 	This parameter has a range of [0, 1].  A value of 0 indicates no
 * 	scattering while a value of 1 indicates no absorption.  Multiple
 * 	scattering tends to dominate materials with high albedos while single
 * 	scattering dominates materials with low albedos.
 *
 * effext:
 * 	Effective extinction defined as follows:
 *
 * 		effext = sqrt(3 * redscatr * absrp)
 *
 * eta:
 * 	Relative index of refraction with the convention that the exterior
 * 	medium's index of refraction is given in the numerator while the
 * 	interior medium's index of refraction is given in the denominator.
 * 	The exterior is defined as the medium that the light is currently in
 * 	while the interior is the medium that the light is entering.
 *
 * ext:
 * 	Extinction coefficient that specifies the probability per unit distance
 * 	that an absorption or scattering event will occur:
 *
 * 		ext = absrp + scatr
 *
 * 	It has the same range as absrp and scatr ([0, +inf]).
 *
 * g:
 * 	Phase function parameter.  This shader uses the Henyey-Greenstein phase
 * 	function.  This parameter has a range of [-1, 1], where -1 corresponds
 * 	to full back scattering, 1 corresponds to full forward scattering, and
 * 	0 corresponds to isotropic scattering.
 *
 * ld:
 * 	Mean-free-path length for light travelling through the surface.
 *
 * rd:
 * 	Diffuse reflectance that relates radiant exitance per unit incident flux
 * 	that will be observed due to a narrow incident beam.  The radiant
 * 	exitance due to multiple scattering is obtained by convolving Rd(r) with
 * 	the irradiance E(x) over a surface:
 *
 * 		Mo(xo) = int { Rd(||xo - xi||) * E(xi) } dA
 *
 * redalbedo:
 * 	Reduced albedo that takes into account the phase parameter:
 *
 * 		redalbedo = redscatr / (redscatr + absrp)
 *
 * redscatr:
 * 	Reduced scattering coefficient that takes into account the phase
 * 	parameter:
 * 	
 * 		redscatr = scatr * (1 - g)
 *
 * 	Since the multiple scattering algorithm assumes isotropic scattering,
 * 	the scattering coefficient must be modified to better estimate the
 * 	amount of reflected light.  Multiplying the scattering coefficient by
 * 	(1 - g) has the effect of increasing scattering for backward scattering
 * 	phase parameters and reducing scattering for forward scattering phase
 * 	parameters.
 *
 * scatr:
 * 	Scattering coefficient that specifies the probability per unit distance
 * 	that a scattering event will occur.  Since the probability is per unit
 * 	distance, this parameter has a range of [0, +inf].  Higher values
 * 	correspond to more scattering.  This parameter is a vector and each
 * 	component corresponds to a different wavelength (R, G, B).
 *
 * td:
 * 	Diffuse transmittance that relates radiant exitance per unit incident
 * 	flux that will be observed due to a narrow incident beam.  The radiant
 * 	exitance due to multiple scattering is obtained by convolving Td(r) with
 * 	the irradiance E(x) over a surface:
 *
 * 		Mo(xo) = int { Td(||xo - xi||) * E(xi) } dA
 */

#ifndef __physicalsss__
#define __physicalsss__

#include <math.h>

#define ERROR_COLOR		{ 1.0, 1.0, 1.0 }
#define MAX_SECANT_ITERATIONS	30
#define PASS_MULTI_GLOBAL_CLR	"pass_multi_global_clr"
#define PC_MODE_GENERATE	0
#define PC_MODE_READ		1
#define PC_MODE_WRITE		2

///
/// Description: Return the diffuse fresnel reflectance for the specified
///		 relative index of refraction.
///
/// Parameters:
///	* eta - relative index of refraction
///
float
getDiffuseFresnelReflectance(float eta)
{
    float fdr;

    if (eta < 1.0)
	fdr =
	    -0.43999 +
	    0.7099 / eta -
	    0.3319 / (eta * eta) +
	    0.0636 / (eta * eta * eta);
    else
	fdr =
	    -1.4399 / (eta * eta) +
	    0.7099 / eta +
	    0.6681 +
	    0.0636 * eta;

    return clamp(fdr, 0.0, 1.0);
}

///
/// Description: Return the diffuse fresnel ratio.  This is a convenience
///		 function that guards against divide-by-zero and provides
///		 clamping.
///
/// Parameters:
///	* eta - relative index of refraction
///
float
getDiffuseFresnelRatio(float eta)
{
    float fdr;

    fdr = clamp(getDiffuseFresnelReflectance(eta), 1e-6, 1.0 - 1e-6);

    return (1.0 + fdr) / (1.0 - fdr);
}

///
/// Description: Returns the evaluation of the Henyey-Greenstein phase function
///		 for the specified directions and phase parameter.
///
/// Parameters:
///	* wi - incoming direction
///	* wo - outgoing direction
///	* g - phase parameter
float
getHenyeyGreensteinPhase(vector wi; vector wo; float g)
{
    return (1.0 - g * g) /
	(2.0 * pow(1.0 - 2.0 * g * dot(wi, wo) + g * g, 1.5));
}

///
/// Description: Computes the integral of the diffuse reflectance from r=0 to
///		 infinity.
///
/// Parameters:
///	* redalbedo - reduced albedo
///	* fdr_ratio - diffuse fresnel ratio
///
float
integrateBrdfRd(float redalbedo; float fdr_ratio)
{
    float a, b, c, d;

    a = sqrt(3.0 * (1.0 - redalbedo));
    b = redalbedo / 2.0;
    c = 1.0 + exp(-(4.0 / 3.0) * fdr_ratio * a);
    d = exp(-a);

    return b * c * d;
}

///
/// Description: Computes the integral of the diffuse reflectance from r=0 to
///		 infinity.
///
/// Parameters:
///	* redalbedo - reduced albedo
///	* fdr_ratio - diffuse fresnel ratio
///
vector
integrateBrdfRdV(vector redalbedo; float fdr_ratio)
{
    vector a, b, c, d;

    a = sqrt(3.0 * ({ 1.0, 1.0, 1.0 } - redalbedo));
    b = redalbedo * 0.5;
    c = { 1.0, 1.0, 1.0 } + exp(-(4.0 / 3.0) * fdr_ratio * a);
    d = exp(-a);

    return b * c * d;
}

vector
vop_sss_opacity(
    vector	x;
    vector	n;
    vector	t;
    vector	absrp;
    vector	scatr)
{
    vector _t, _absrp, _scatr, nmax;
    float dist;
    int hit;

    if (dot(n, t) > 0.0)
	return 0.0;

    // Clamp/normalize input.
    _t = normalize(t);
    _absrp = max({ 0.0, 0.0, 0.0 }, absrp);
    _scatr = max({ 0.0, 0.0, 0.0 }, scatr);

    // Trace through the surface until we hit the back side.
    hit = trace(
	x, _t, Time, "pipeline", "displacement", "scope", "scope:self",
	"SID", SID, "N", nmax, "ray:length", dist);

    if (!hit || dot(nmax, _t) < 0.0)
	return 0.0;

    // Determine the opacity of the attenuated light.
    return { 1.0, 1.0, 1.0 } - exp(-(_absrp + _scatr) * dist);
}

///
/// Description: Returns the reduced albedo for the specified diffuse
///		 reflectance and relative index of refraction.  threshold is
///		 used as a stopping criteria for a secant root finding
///		 algorithm.
///
/// Parameters:
///	* rd - diffuse reflectivity
///	* eta - relative index of refraction
///	* threshold - a small number used as an error threshold
///
vector
vop_sss_reduced_albedo(const vector diff;
		       const float eta;
		       const float threshold)
{
    // Change to reflectance and clamp to avoid floating point precision
    // errors in the resulting reduced albedo.
    vector rd = clamp(diff * diff, 0, 0.99);

    // Clamp/normalize input.
    float _eta = max(0.0, eta);

    float fdr_ratio = getDiffuseFresnelRatio(1.0 / _eta);

    vector x;
    for (int icomp = 0; icomp < 3; ++icomp)
    {
	float x1 = 0.0;
	float x2 = 1.0;
	float f1 = integrateBrdfRd(x1, fdr_ratio) - rd[icomp];
	float f2 = integrateBrdfRd(x2, fdr_ratio) - rd[icomp];
	int niterations = 0;

	while (niterations < MAX_SECANT_ITERATIONS &&
	   (f2 > threshold || f2 < -threshold) && (f2 - f1 != 0.0))
	{
	    float x3 = x2 - f2 * ((x2 - x1) / (f2 - f1));

	    x1 = x2;
	    x2 = clamp(x3, 0.0, 1.0);
	    f1 = f2;
	    f2 = integrateBrdfRd(x2, fdr_ratio) - rd[icomp];
	    niterations++;
	}

	x[icomp] = x2;
    }

    return x;
}

///
/// Description: Returns the diffuse reflectivity for the specified reduced
///		 albedo and relative index of refraction.
///
/// Parameters:
///	* redalbedo - reduced albedo
///	* eta - relative index of refraction
///
vector
vop_sss_diffuse_reflectivity(const vector redalbedo; const float eta)
{
    // Clamp/normalize input.
    float _eta = max(0.0, eta);

    float fdr_ratio = getDiffuseFresnelRatio(1.0 / _eta);

    // Compute the reflectance
    vector rd = (redalbedo/2) *
	(1 + exp(-(4/3) * fdr_ratio * sqrt(3*(1-redalbedo)))) *
	exp(-sqrt(3*(1-redalbedo)));

    // Change to reflectivity
    return sqrt(rd);
}

#endif

