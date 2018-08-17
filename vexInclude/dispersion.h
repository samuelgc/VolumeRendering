#ifndef __dispersion_utils_h__
#define __dispersion_utils_h__

#include <cie_cmf.h>

#define FRAUN_D 586.17
#define DISPERSION_PIVOT 0.234812930625 // pow(1.0 - index_fraction(FRAUN_D),2.0)
					// centers abbe pivot about Fraunhofer D line

float index_fraction( const float wavelength ) {
    return (wavelength - CMF_WLSTART) / (CMF_WLEND - CMF_WLSTART);
}

float
disperse_ior(const float ior, dispersion, wavelength)
{
	float delta_ior = 1.0f - index_fraction(wavelength);
	delta_ior *= delta_ior;
	return ior + (delta_ior - DISPERSION_PIVOT) * dispersion;
}

float
cmf_importance( float pdf;
                float in_sx
                )
{
    // slightly fudged to allow more samples around low-brightness regions
    int tableSize = 32;
    float ltable[] = {
        0.0916, 0.1036, 0.1814, 0.3660, 0.6225, 0.7225, 0.7547, 0.7856,
        0.8636, 1.0878, 1.4875, 1.9784, 2.5091, 2.7602, 2.9183, 2.9147,
        2.7677, 2.5146, 2.0991, 1.5986, 0.9479, 0.5736, 0.3343, 0.2049,
        0.1400, 0.1114, 0.0993, 0.0943, 0.0920, 0.0920, 0.0920, 0.0912
     };

    float r0 = in_sx * tableSize;
    int sel = 0;
    for (sel= 0; sel< tableSize; ++sel)
    {
        r0 -= ltable[sel];
        if (r0 < 0) break;
    }
    pdf = ltable[sel];
    return (sel - r0/ltable[sel]) / tableSize;
}

#endif
