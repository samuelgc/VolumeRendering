#ifndef __bsdf_h__
#define __bsdf_h__

#include "voplib.h"
#include "pbr.h"
#include "ggx_utils.h"

// To avoid numerical issues due to division by 0
#define SPEC_MIN_WIDTH 1e-5

// To avoid numerical issues due taking logarithms of a value close to 0.
// We'll blend with an exponent of 1 when the cosine is less than this
// value.
#define SPEC_MIN_COSINE2 0.1
#define LOG4 1.3862943611198906

// Convert a highlight angle to an exponent to be passed to the phonglobe()
// operation.
float
pbrspecular_toexponent(float width)
{
    float       costheta2;
    float       costheta2_clamped;
    float       exponent;

    costheta2 = max(width, SPEC_MIN_WIDTH);
    costheta2 = 1-costheta2*costheta2;
    costheta2_clamped = max(costheta2, SPEC_MIN_COSINE2);

    exponent = -2*log(2)/log(costheta2_clamped) + 1;
    if (costheta2 < 0.1)
    {
        exponent = lerp(1.0, exponent, costheta2 / SPEC_MIN_COSINE2);
    }

    return exponent;
}

float
pbrspecular_rough_to_exponent(float roughness)
{
    return (2.0/(roughness*roughness)-2.0);
}

float
pbrspecular_rough_to_angle(float roughness)
{
    float invrough = 2.0/roughness;

    float costheta = exp(LOG4/(-invrough*invrough+6));

    //at 0.810805 roughness we have 90 degrees.
    //guard against higher values here. should maybe exit earlier.
    return costheta > 1 ? 0.5 * M_PI : acos(costheta);
}

float eta_to_fzero(float eta)
{
    return (eta*eta - 2.0*eta + 1.0) / (eta*eta + 2.0*eta + 1.0);
}

vector combined_fresnel_reflect(float eta, kr; vector F0)
{
    // This is to be used with fresnel blending.
    //
    // F0 is repurposed F0 as a tint color.
    //
    // For refraction, we simply multiply with F0.
    //
    // For reflection, at facing angles we use F0 multiplied by the
    // facing reflectivity corresponding eta value. At grazing angles
    // we use full, untinted reflectivity.

    float fresf0 = eta_to_fzero(eta);

    return lerp(fresf0*F0, {1,1,1}, fit(kr, fresf0, 1.0, 0.0, 1.0));
}

bsdf get_bsdf(string model; vector nNg, nN, nI, utan, vtan, F0;
	      float eta, rough, aniso, aniso_dir;
	      int masking, thinfilm, fresblend;
	      float reflect, refract, dispersion;
	      export float handled_energy; string reflectlabel, refractlabel)
{
    // combined, reflection and refraction BSDFs
    bsdf f, fr, ft = bsdf();
    handled_energy = 0.0;

    vector utan_ = utan;    // temporary tangent vectors
    vector vtan_ = vtan;

    float kr;		    // reflection intensity
    float kt;		    // transmission intensity
    vector nR;		    // normalized reflection vector
    vector nT;		    // normalized transmission vector

    // check whether we have perfect mirror specular
    int perfect_spec = rough <= SPEC_MIN_WIDTH;

    if (model != "ggx")// || perfect_spec)
    {
	if(model != "ggx" && thinfilm)
	{
	    thinfresnel(nI, nN, eta, kr, kt, nR, nT);
	}
	else
	{
	    fresnel(nI, nN, eta, kr, kt, nR, nT);
	}
    }

    float dotNI = dot(nN, nI);

    vector flippedNR = select(dotNI > 0, -nN, nN);
    vector flippedNT = select(dotNI < 0, -nN, nN);

    if (!perfect_spec || model == "ggx")
    {
        float       uexp, vexp;
        float       rough2;

        if(model == "blinn" || model == "ggx")
        {
            // square roughness to get a more linear mapping
            rough2 = rough*rough;
        }
        else
        {
            rough2 = rough;
        }

        uexp = vexp = pbrspecular_rough_to_exponent(rough2);

        if ((model != "cone" && abs(aniso) > SPEC_MIN_WIDTH) && !perfect_spec)
        {
            if (aniso_dir != 0)
            {
                matrix3 rot = ident();
                rotate(rot, M_PI * aniso_dir, cross(utan_, vtan_));
                utan_ *= rot;
                vtan_ *= rot;
            }



            if (aniso < 0)
                uexp *= pbrspecular_toexponent(1+(aniso));
            else
                vexp *= pbrspecular_toexponent(1-(aniso));

            if (model == "blinn")
            {
                fr = ashikhmin(flippedNR, uexp, vexp, utan_, vtan_, "label", reflectlabel);
            }
            else if (model == "ggx")
            {
                // make orthonormal
                vtan_ = normalize(cross(nN,utan_));
                utan_ = normalize(cross(nN,vtan_));

                int flip = aniso < 0;

                float aspect = sqrt(1.0-abs(aniso)*0.9);
                float alphax = max(0.001, rough2/aspect);
                float alphay = max(0.001, rough2*aspect);

		f = cvex_bsdf(
			       "oplib:/Shop/ggx_eval?Shop/ggx_eval",
			       "oplib:/Shop/ggx_sample?Shop/ggx_sample",
			       "label", reflectlabel + " " + refractlabel,
			       "ng", nNg,
			       "nn", nN,
			       "xg", select(flip, vtan_, utan_),
			       "yg", select(flip, utan_, vtan_),
			       "F0", F0,
			       "alphax", alphax,
			       "alphay", alphay,
			       "masking", masking,
			       "fresblend", fresblend,
			       "eta", eta,
			       "reflect", reflect,
			       "refract", refract,
			       "reflectmask", bouncemask(reflectlabel),
			       "refractmask", bouncemask(refractlabel),
			       "dispersion", dispersion);
            }
            else
            {
                fr = phonglobe(flippedNR, nR, uexp, vexp, utan_, vtan_, "label", reflectlabel);
                ft = phonglobe(flippedNT, nT, uexp, vexp, utan_, vtan_, "label", refractlabel);
            }
        }
        else
        {
            if (model == "cone")
            {
                float newangle = pbrspecular_rough_to_angle(rough2);
                fr = cone(flippedNR, nR, newangle, "label", reflectlabel);
                ft = cone(flippedNT, nT, newangle, "label", refractlabel);
            }
            else if (model == "blinn")
            {
                fr = blinn(flippedNR, uexp, "label", reflectlabel);
            }
            else if(model == "ggx")
            {
		// create arbitrary tangent vectors for the isotropic case.
		matrix3 to_world = dihedral({0.0,0.0,1.0}, nN);
		utan_ = {1,0,0} * to_world;
		vtan_ = {0,1,0} * to_world;

                int flip = aniso < 0;

                float alpha = max(0.001, rough2);

		f = cvex_bsdf(
			       "oplib:/Shop/ggx_eval?Shop/ggx_eval",
			       "oplib:/Shop/ggx_sample?Shop/ggx_sample",
			       "label", reflectlabel + " " + refractlabel,
			       "ng", nNg,
			       "nn", nN,
			       "xg", select(flip, vtan_, utan_),
			       "yg", select(flip, utan_, vtan_),
			       "F0", F0,
			       "alphax", alpha,
			       "alphay", alpha,
			       "masking", masking,
			       "fresblend", fresblend,
			       "eta", eta,
			       "reflect", reflect,
			       "refract", refract,
			       "reflectmask", bouncemask(reflectlabel),
			       "refractmask", bouncemask(refractlabel),
			       "dispersion", dispersion);
            }
            else
            {
                fr = phonglobe(flippedNR, nR, uexp, "label", reflectlabel);
                ft = phonglobe(flippedNT, nT, uexp, "label", refractlabel);
            }
        }
	if(model != "ggx")
	{
	    f  *= 1.0 / luminance(albedo(f));
	    fr *= 1.0 / luminance(albedo(fr));
	    ft *= 1.0 / luminance(albedo(ft));
	}
    }
    else // we have a tiny roughness value, use the efficient specular BSDF
    {
	fr = specular(nR, "label", reflectlabel);
	ft = specular(nT, "label", refractlabel);
    }

    // Fresnel and/or Schlick-based blending.
    //
    // Microfacet models do fresnel/schlick blending internally, so we only
    // need to do it here if we have prefect specular/no roughness, in which
    // case we used the 'specular' BSDF instead.
    //
    // For non-microfacet BSDFs, always do the blending at this stage.
    if (model != "ggx")// || perfect_spec)
    {
	vector trans_int;
	vector refl_int;
	int handle_both = reflect > 0 && refract > 0;
	if(fresblend || handle_both)
	{
	    if(fresblend)
	    {
		kt = 1.0 - (reflect*kr);
		handled_energy = lerp(reflect * kr, 1, refract);
	    }
	    else
	    {
		kr = 1;
		kt = 1;
		handled_energy = refract;
	    }

	    refl_int = reflect * combined_fresnel_reflect(eta, kr, F0);
	    trans_int = refract * kt;

	    f = trans_int * ft + refl_int * fr;

	}
	else
	{
	    float schlick = schlick_f(abs(dot(nN, nI)));
	    trans_int = refract;
	    refl_int = reflect * F0 + (1.0-F0) * schlick;
	    f = select(refract > 0, ft * trans_int, fr * refl_int);
	    handled_energy = refract;
	}
    }
    else
    {
	float est_eta = select(fresblend, 1.0/eta, 0.0);
	vector estimate = ggx_albedo_estimate(dot(nN, nI), rough, est_eta, F0);
	handled_energy = lerp(reflect * luminance(estimate), 1.0, refract);
    }

    return f;
}

bsdf get_bsdf(string model, label; vector nNg, nN, nI, utan, vtan, F0; float eta, rough, aniso, aniso_dir; int masking, thinfilm, fresblend)
{
    float handled_energy;
    float refract = (bouncemask(label) & PBR_REFRACT_MASK) > 0;
    float reflect = refract == 0;
    float dispersion = 0;
    return get_bsdf(model, nNg, nN, nI, utan, vtan, F0, eta, rough, aniso, aniso_dir, masking, thinfilm, fresblend, reflect, refract, dispersion, handled_energy, label, label);
}

bsdf get_bsdf(string model, label; vector nNg, nN, nI, utan, vtan; vector F0; float eta, rough, aniso, aniso_dir; int masking, thinfilm)
{
    return get_bsdf(model, label, nNg, nN, nI, utan, vtan, F0, eta, rough, aniso, aniso_dir, masking, thinfilm, 0 /*fresblend*/);
}

bsdf get_bsdf(string model, label; vector nNg, nN, nI, utan, vtan; vector F0; float eta, rough, aniso, aniso_dir; int masking)
{
    return get_bsdf(model, label, nNg, nN, nI, utan, vtan, F0, eta, rough, aniso, aniso_dir, masking, 0 /*thinfilm*/, 0 /*fresblend*/);
}

#endif
