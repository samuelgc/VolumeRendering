//
// From course notes for "Advanced RenderMan 3" Siggraph 2001 by Matt Pharr
// 
// Surface shader that implements a shading model that should have a visual 
// appearence generall similar to that of skin.  Based on phenomenological 
// information about skin reflectance from Hanrahan and Krueger, 
// "Reflection from layered surfaces due to subsurface scattering", 
// proceedings of Siggraph 1993.
//

/* Evaluate the Henyey-Greenstein phase function for two vectors with
   an asymmetry value g.  v1 and v2 should be normalized and g should 
   be in the range (-1, 1).  Negative values of g correspond to more
   back-scattering and positive values correspond to more forward scattering.
*/
float
phase(vector v1, v2; float g)
{
    float costheta = dot(-v1, v2);
	float g2 = g*g;
    return (1.0 - g2) / pow(1.0 + g2 - 2.*g*costheta, 1.5);
}

/* Compute a the single-scattering approximation to scattering from
   a one-dimensional volumetric surface.  Given incident and outgoing
   directions wi and wo, surface normal n, asymmetry value g (see above),
   scattering albedo (between 0 and 1 for physically-valid volumes),
   and the depth of the volume, use the closed-form single-scattering
   equation to approximate overall scattering.
*/
float
singleScatter(vector in, out, nn; float g, albedo, depth) 
{
    float win = abs(dot(in, nn));
    float won = abs(dot(out, nn));
    float offset = phase(out, in, g);

    float scatter = albedo * offset/(win + won);
    scatter *= 1.0 - exp(-(1.0/win + 1.0/won) * depth);

    return scatter;
}

vector
efresnel(vector ii, nn; float eta; float Kr, Kt;) 
{
    vector R, T;
    fresnel(ii, nn, eta, Kr, Kt, R, T);
    Kr = smooth(0.0, 0.5, Kr);
    Kt = 1.0 - Kr;
    return normalize(T);
}
