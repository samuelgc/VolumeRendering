/*
* PROPRIETARY INFORMATION.  This software is proprietary to
* Side Effects Software Inc., and is not to be reproduced,
* transmitted, or disclosed in any way without written permission.
*
* Produced by:
*       Side Effects Software Inc
*       477 Richmond Street West
*       Toronto, Ontario
*       Canada   M5V 3E7
*       416-504-9876
*
* NAME: ocean.h ( Ocean support functions, VEX )
*
* COMMENTS:     Ocean support functions.
*/
#ifndef __ocean_h_
#define __ocean_h_

#include <math.h>
#include <complex.h>

#define INPUT( num ) (sprintf("opinput:%d", (num) ))

/// Random values in a Guassian/Normal distribution using the Box-Muller
/// transform.
float
grandom(const float u0, u1, sigma)
{
    // Given two independent random variables u0 and u1, use the Box-Muller
    // transformation to generate a Gaussian random distribution with the given
    // standard deviation.
    return cos(u0 * (M_PI*2)) * sqrt(-2 * log(u1)) * sigma;
}

// Flatten 3d vector components to 2d along the given plane.
vector2
twod(const int x, y, z, plane)
{
    if (plane == 0)
        return set(x, y);
    if (plane == 1)
        return set(y, z);
    return set(x, z);
}

// Flatten a 3d vector to 2d along the given plane.
vector2
twod(const vector v; int plane)
{
    if (plane == 0)
        return set(v.x, v.y);
    if (plane == 1)
        return set(v.y, v.z);
    return set(v.x, v.z);
}

// Expand 2d vector components given the original plane.
vector
threed(const int x, y, plane)
{
    if (plane == 0)
        return set(x, y, 0);
    if (plane == 1)
        return set(0, x, y);
    return set(x, 0, y);
}

// Expand a 2d vector to 3d given the original plane.
vector
threed(const vector2 v; const int plane)
{
    if (plane == 0)
        return set(v.x, v.y, 0);
    if (plane == 1)
        return set(0, v.x, v.y);
    return set(v.x, 0, v.y);
}

// Return the value along the planar dimension.
float
planar(const vector a; const int plane)
{
    return getcomp(a, (plane + 2) % 3);
}

// Sets the value along the planar dimension.
void
setplanar(vector a; float value; const int plane)
{
    setcomp(a, value, (plane + 2) % 3);
}

// Returns which plane of the supplied resolution is 2d.
int
plane2d(const vector res)
{
    if (res.z == 1)
        return 0;
    if (res.x == 1)
        return 1;
    // Fallback to assuming non-square axis represents planar.
    if (res.x == res.y)
        return 0;
    if (res.y == res.z)
        return 1;
    return 2;
}

// Returns which plane of the supplied volume is 2d.
int
plane2d(const string opname; const string volname)
{
    vector res = volumeres(opname, volname);
    return plane2d(res);
}

// Returns which plane of the supplied volume is 2d.
int
plane2d(const int geo; const string volname)
{
    return plane2d(INPUT(geo), volname);
}

// Returns which plane is 2d by looking at known ocean volume names.
int
oceanPlane2d(const string opname)
{
    string volumes[] = { "amplitude", "restdisplace.x", "restvel.x", "restderivatives.x" };
    foreach(string vol; volumes)
        if (nametoprim(opname, vol) >= 0)
            return plane2d(opname, vol);
    // Default to xz.
    return 2;
}

// Expand the 2d derivatives in dx to a full 3x3 Jacobian,
// given the 2d plane.
matrix3
jacobian3d(const vector dx; const int plane)
{
    if (plane == 0)
        return set( dx.x,   dx.y,   0,
                    dx.y,   dx.z,   0,
                    0,      0,      0);
    if (plane == 1)
        return set( 0,      0,      0,
                    0,      dx.x,   dx.y,
                    0,      dx.y,   dx.z);
    return set( dx.x,   0,  dx.y,
                0,      0,  0,
                dx.y,   0,  dx.z);
}

// Collapse the 3d Jacobian to elements of 
// 2d symmetric sub-matrix, given the 2d plane.
// (i.e. the converse of above function)
vector
jacobian2d(const matrix3 J; const int plane)
{
    if (plane == 0)
        return set(J.xx, J.xy, J.yy);
    if (plane == 1)
        return set(J.yy, J.yz, J.zz);
    return set(J.xx, J.xz, J.zz);
}

// Convert from linear point / voxel number to 3-d voxel index, assuming
// volume is 2d.
vector
lineartoidx2d(const int geo; const string volname; const int n)
{
    vector res = volumeres(geo, volname);
    int plane = plane2d(res);
    vector r2 = twod(res, plane);
    int num = n;
    int x = num % (int)r2.x;
    num -= x;
    num /= (int)r2.x;
    int y = num % (int)r2.y;
    return threed(x, y, plane);
}

// Compute the 2d wavevector for the given voxel index and resolution.
// If normalized is true, ignore gridsize and return a value based only on
// index and res, suitable for random number generation.
vector2
wavevector(const vector2 voxidx, res;
           const float gridsize; const int normalized)
{
    vector2 k =  (voxidx - floor(res / 2));
    if (!normalized)
        k *= M_TWO_PI / gridsize;
    return k;
}

// Returns the smallest wavelength represented in the given spectrum.
float minwavelength(const vector2 res; const float gridsize; const int downsampled)
{
    // Use the zero'th coefficient, i.e. largest negative wavevector.
    float l = M_TWO_PI / length(wavevector(set(0, 0), res, gridsize, 0));
    // Each level of downsampling doubles the min wavelength.
    if (downsampled)
        l *= pow(2, downsampled);
    return l;
}

// Returns true if the current idx is for the zero wavevector.
int
iswavezero(const vector2 idx, res)
{
    return (idx == res/2);
}

// Compute the amplitude for the given wavevector from the Phillips statistical
// model.
float
phillipsspectrum(const vector2 wavevector;
                 const float winddir, windspeed, gravity,
                 scale, windbias, windmove, mink)
{
    float k = length(wavevector);
    float L = windspeed * windspeed / gravity;
    float k2 = k * k;
    float k2L2 = k2 * L * L;
    float k2l2 = k2 * mink * mink;
    float k4 = k2 * k2;
    float P = scale * exp(-1 / k2L2 - k2l2) / k4;
    if (windbias > 0)
    {
        // Compute wind power as amount wave aligned with wind.
        vector2 windvector = set(cos(winddir), sin(winddir));
        float windp = dot(windvector, wavevector) / k;
        // Wave moving in opposite direction from wind, reduce wind power.
        if (windp < 0 && windmove > 0)
            windp *= clamp((1 - windmove), 0, 1);
        // Increase amplitudes for waves aligned with wind in either direction.
        P *= pow(abs(windp), windbias);
    }
    return sqrt(P) * M_SQRT1_2;
}

// Identify the 2d plane of the grid, project the voxel
// index and resolution into 2d vectors and calculate gridsize.
void
project2dplane(const int ix, iy, iz, resx, resy, resz;
               const vector dPdx, dPdy, dPdz;
               vector2 idx, res; float gridsize; int plane)
{
    plane = plane2d(set(resx, resy, resz));
    idx = twod(ix, iy, iz, plane);
    res = twod(resx, resy, resz, plane);
    gridsize = max(dPdx.x * resx, dPdy.y * resy, dPdz.z * resz);
}

// Clamp the input vector to zero and the bottom of the last voxel of the input volume
// along the specified plane.  This allows us to sample from a volume with repeat
// boundaries along the horizontal dimensions without repeating along the depth.
vector
clampToZeroMaxDepth(const int geo; const string volname; const vector p; const int plane)
{
    vector minb, maxb, res;
    getbbox(geo, minb, maxb);
    res = volumeres(geo, volname);
    vector size = maxb - minb;
    // Bottom of last voxel.
    float maxdepth = planar(minb + size / (2 * res), plane);
    maxdepth = max(maxdepth, planar(p, plane));
    vector clampedp = p;
    setplanar(clampedp, min(0, maxdepth), plane);
    return clampedp;
}

// Compute the frequency and phase for the given wavevector based on
// graviy, depth, and the dispersion relation.
void
wavefrequencyphase(const vector2 voxidx, res;
                   const float gridsize, gravity;
                   const int seed;
                   const float depth;
                   const int loop;
                   const float loopperiod;
                   float frequency, phase)
{
    // Calc wavevector and wavenumber.
    vector2 k = wavevector(voxidx, res, gridsize, 0);
    float kl = length(k);
    // Calc wave frequency with dispersion.
    frequency = sqrt(tanh(depth * kl) * kl * gravity);

    if (loop)
    {
        // Quantize frequencies to repeat at loopperiod.
        float quantperiod = M_TWO_PI / max(0.01, loopperiod);
        frequency = int(frequency / quantperiod) * quantperiod;
    }

    // Base the phase random number on the normalized wave vector
    // so that we always get same result for a given wavevector as we
    // change gridsize and/or resolutions.
    vector2 kn = wavevector(voxidx, res, gridsize, 1);
    // Calc phase as random number from 0..2pi.
    phase = fit01(random(set(kn.x, kn.y, seed)), 0, M_TWO_PI);
}

// Compute the amplitude for the given wavevector based on
// wind parameters fed into the Phillips spectrum model, with random scaling
// by the given distribution. (1 = uniform, 2 = Gaussian, 3 = Log Normal)
float
waveamplitude(const vector2 voxidx, res;
              const float gridsize, scale, gravity;
              const int seed;
              const float windspeed, dirmove, dirbias, winddir, mink;
              const int distribution)
{
    // Use the Phillips spectrum model for amplitudes.
    vector2 k = wavevector(voxidx, res, gridsize, 0);
    float amp = phillipsspectrum(k, winddir, windspeed, gravity,
                                scale, dirbias, dirmove, mink);

    // Base the amplitude random scales on the normalized wave vector
    // so that we always get same result for a given wavelength as we
    // change resolutions.
    if (distribution > 0)
    {
        vector2 kn = wavevector(voxidx, res, gridsize, 1);
        float r1 = random(set(kn.x, seed, kn.y));
        float r2 = random(set(kn.x, seed * 123, kn.y));

        if (distribution == 1) // Uniform
            amp *= fit01(r1, -1, 1);
        else if (distribution == 2) // Gaussian
            amp *= grandom(r1, r2, 1);
        else if (distribution == 3) // Log Normal
            amp *= exp(grandom(r1, r2, 1));
    }

    // Ensure we have no DC offset from wave zero.
    if (iswavezero(voxidx, res))
        amp = 0;
    return amp;
}

// Compute a 0-1 value to remap the wave amplitude based on the current
// wavenumber relative to the max wavenumber.
float
remapamplitude(const vector2 voxidx, res; const float gridsize;
               const int maxremap)
{
    int maxres = (int)pow(2, maxremap);
    vector2 k = wavevector(voxidx, res, gridsize, 0);
    vector2 maxk = wavevector(0, maxres, gridsize, 0);
    return length(k) / length(maxk);
}

// Returns the offset into a hi-res volume when downsampling.
vector2
downsampleOffset(const vector2 hires; int downsample)
{
    vector2 offset = hires - hires / pow(2, downsample);
    return offset / 2;
}


// Calculate the positive and negative components of the Fourier coefficient
// for Y displacement:
// Y = 1/2 * exp(i * (phase_i_j - frequency * time)) * amplitude_i_j +
//     1/2 * exp(-i * (phase_-i_-j - frequency * time)) * amplitude_-i_-j
void
displacementY(const int geo; const string ampname, phasename;
              const vector2 idx, res; int plane;
              const float amplitude, frequency, phase, time;
              complex Ypos, Yneg)
{
    // Positive component.
    Ypos =  cmultreal(cexpimag(phase - frequency * time), 0.5 * amplitude);

    // Look up values comprising the negative component.
    vector negfreqidx = threed(res - idx, plane);
    float ampneg = volumeindex(geo, ampname, negfreqidx);
    float phaseneg = volumeindex(geo, phasename, negfreqidx);
    // Frequency is function of wavenumber so is same
    // for positive and negative wavevectors.
    float freqneg = frequency;
    // Negative component.
    Yneg = cmultreal(cexpimag(-(phaseneg - freqneg * time)), 0.5 * ampneg);
}

// Calculate the Fourier coefficient for dY / dt from the positive and negative
// components of Y:
// V = dY / dt = i * -frequency * Ypos + i * frequency * Yneg
complex
velocityY(const complex Ypos, Yneg; const float frequency)
{
    return cadd(cmultimag(Ypos, -frequency), cmultimag(Yneg, frequency));
}

// Calculate the Fourier coefficient for horizontal displacement
// from the vertical Y displacement and a tunable chop parameter.
// NOTE: this gives horizontal velocity from vertical velocity as well.
// (X, Z) = Y * (i * chop * normalize(k))
void
displaceXZ(const vector2 idx, res; const float gridsize, chop;
           const complex Y; complex X, Z)
{
    if (chop > 0)
    {
        vector2 k = wavevector(idx, res, gridsize, 0);
        k = chop * normalize(k);
        X = cmultimag(Y, k.x);
        Z = cmultimag(Y, k.y);
    }
    else
    {
        X = Z = CMPLX(0, 0);
    }
}

// Calculate the Fourier coefficients for the derivatives of the
// X and Z displacement (where D is the horizontal displacement)
// (dDx / dx, dDx / dz) = X * (i * k)
// (dDz / dx, dDz / dz) = Z * (i * k)
void
derivativesXZ(const vector2 idx, res; const float gridsize;
           const complex X, Z; complex dxx, dxz, dzz)
{
    vector2 k = wavevector(idx, res, gridsize, 0);
    dxx = cmultimag(X, k.x);
    dxz = cmultimag(X, k.y);
    dzz = cmultimag(Z, k.y);
}

// Compute cusp and cusp direction values given the XZ derivatives in dx.
void
computeCusp(const vector dx; float cusp; vector cuspdir; const int plane)
{
    // Compute min eigenvalues and eigenvectors of Jacobian matrix.
    // Full jacobian of horizontal displacement is
    // J([x,y,z]) + J = I + J 
    float Jxx = dx.x + 1;   // Jxx = dXdx + 1
    float Jxz = dx.y;       // Jxz = dXdz (== dZdx)
    float Jzz = dx.z + 1;   // Jzz = dZdx + 1

    // Calculate the minimum eigenvalue and eigenvector of the symmetrix 2x2
    // Jacobian matrix given the real spatial partial derivatives.
    // J = 1/2 * [(Jxx + Jzz) - sqrt((Jxx - Jzz)^2 + 4 * Jxz^2))]
    // q = normalize(1, (J - Jxx) / Jxz)
    float mineig;
    vector2 mineigvec;
    float K = Jxx - Jzz;
    K *= K;
    K += 4 * Jxz * Jxz;
    mineig = 0.5 * ((Jxx + Jzz) - sqrt(K));
    cusp = 1 - mineig;
    // If mineig is almost 1, too flat to get useful cusp / cuspdir.
    if (abs(cusp) < 1e-5)
    {
        cusp = 0;
        cuspdir = 0;
        return;
    }    

    mineigvec = normalize(set(1, (mineig - Jxx) / Jxz));
    cuspdir = threed(mineigvec, plane);
}

// Compute cusp and cusp direction values given the 3x3 Jacobian.
void
computeCusp(const matrix3 J; float cusp; vector cuspdir; const int plane)
{
    vector offplane = 0;
    setplanar(offplane, 1, plane);
    // If J has a zero offplane column, we can calc faster in 2d.
    if (length2(offplane * J) < 1e-4)    
    {
        computeCusp(jacobian2d(J, plane), cusp, cuspdir, plane);
        return;
    }

    // Off-axis rotations involved, so we need full 3d calc.
    // Full jacobian of horizontal displacement is
    // J([x,y,z]) + J = I + J 
    matrix3 Im = 1;
    matrix3 Jx = Im + J;
    int nroots;
    vector real, imag;
    eigenvalues(nroots, Jx, real, imag);
    // If nroots is 1, Jx is close enough to I
    // to bail out.
    if (nroots == 1)
    {
        cusp = 0;
        cuspdir = 0;
        return;
    }
    // With three roots, mineig is always first.
    // (We should always have three since J *should be*
    // symmetric.)
    float mineig = real.x;
    cusp = 1 - mineig;
    // If mineig is almost 1, too flat to get useful cusp / cuspdir.
    if (abs(cusp) < 1e-5)
    {
        cusp = 0;
        cuspdir = 0;
        return;
    }
    // Use Cayley-Hamilton theorem with other eigenvalues
    // to compute matrix whose columns are the eigenvector for
    // mineig.
    matrix3 e = (Jx - real.y * Im) * (Jx - real.z * Im);
    cuspdir = normalize(set(e.xx, e.yx, e.zx));
}

int
primByName(const string opname; const string name; const int which)
{
    return findattribval(opname, "primitive", "name", name, which);
}

int
primByName(const int geo; const string name; const int which)
{
    return primByName(INPUT(geo), name, which);
}

int
nprimsByName(const string opname; const string name)
{
    return findattribvalcount(opname, "primitive", "name", name);
}

int
nprimsByName(const int geo; const string name)
{
    return nprimsByName(INPUT(geo), name);
}

float
accumulatevolumes(const string opname; const string vol; const vector p; const int start)
{
    int nvol = nprimsByName(opname, vol);
    float result = 0;
    for (int i=start; i < nvol; i++)
    {
        int prim;
        prim = primByName(opname, vol, i);
        if (prim >= 0)
            result += volumesample(opname, prim, p);
    }
    return result;
}

// Accumulate float values from the named volume, beginning from start.
float
accumulatevolumes(const int geo; const string vol; const vector p; const int start)
{
    return accumulatevolumes(INPUT(geo), vol, p, start);
}

// Accumulate float values from the named volume, beginning from zero.
float
accumulatevolumes(const int geo; const string vol; const vector p)
{
    return accumulatevolumes(INPUT(geo), vol, p, 0);
}

// Accumulate vector values from the named vector volumes, beginning from start.
vector
accumulatevolumesv(const string opname; const string x, y, z; const vector p; const int start)
{
    int nvol = nprimsByName(opname, x);
    vector result = 0;
    for (int i=start; i < nvol; i++)
    {
        int prim;
        prim = primByName(opname, x, i);
        if (prim >= 0)
            result.x += volumesample(opname, prim, p);
        prim = primByName(opname, y, i);
        if (prim >= 0)
            result.y += volumesample(opname, prim, p);
        prim = primByName(opname, z, i);
        if (prim >= 0)
            result.z += volumesample(opname, prim, p);
    }
    return result;
}

// Accumulate vector values from the named vector volumes, beginning from start.
vector
accumulatevolumesv(const int geo; const string x, y, z; const vector p; const int start)
{
    return accumulatevolumesv(INPUT(geo), x, y, z, p, start);
}

// Accumulate vector values from the named vector volumes, beginning from zero.
vector
accumulatevolumesv(const int geo; const string x, y, z; const vector p)
{
    return accumulatevolumesv(INPUT(geo), x, y, z, p, 0);
}

// Accumulate vector values from the vector volumes that start with prefix,
// beginning from start.
vector
accumulatevolumesv(const string opname; const string prefix; const vector p; int start)
{
    string x, y, z;
    x = sprintf("%s.x", prefix);
    y = sprintf("%s.y", prefix);
    z = sprintf("%s.z", prefix);
    return accumulatevolumesv(opname, x, y, z, p, start);
}

// Accumulate vector values from the vector volumes that start with prefix,
// beginning from zero.vector
vector
accumulatevolumesv(const string opname; const string prefix; const vector p)
{
    return accumulatevolumesv(opname, prefix, p, 0);
}

// Accumulate vector values from the vector volumes that start with prefix,
// beginning from start.
vector
accumulatevolumesv(const int geo; const string prefix; const vector p; int start)
{
    return accumulatevolumesv(INPUT(geo), prefix, p, start);
}

// Accumulate vector values from the vector volumes that start with prefix,
// beginning from zero.
vector
accumulatevolumesv(const int geo; const string prefix; const vector p)
{
    return accumulatevolumesv(INPUT(geo), prefix, p, 0);
}

// Version of volumeindex that allows selecting which scalar volume to index.
float
volumeindex(const string opname; const string vol; const vector idx; int which)
{
    float result = 0;
    int prim;
    prim = primByName(opname, vol, which);
    if (prim >= 0)
        result = volumeindex(opname, prim, idx);
    return result;
}

// Version of volumeindexv that allows selecting which vector volume to index.
vector
volumeindexv(const string opname; const string prefix; const vector idx; int which)
{
    string x, y, z;
    x = sprintf("%s.x", prefix);
    y = sprintf("%s.y", prefix);
    z = sprintf("%s.z", prefix);
    vector result = 0;
    int prim;
    prim = primByName(opname, x, which);
    if (prim >= 0)
        result.x = volumeindex(opname, prim, idx);
    prim = primByName(opname, y, which);
    if (prim >= 0)
        result.y = volumeindex(opname, prim, idx);
    prim = primByName(opname, z, which);
    if (prim >= 0)
        result.z = volumeindex(opname, prim, idx);
    return result;
}

vector
volumeindexv(const int geo; const string prefix; const vector idx; int which)
{
    return volumeindexv(INPUT(geo), prefix, idx, which);
}

// Version of volumesample that allows selecting which volume to index.
float
volumesample(const int geo; const string vol; const vector p; int which)
{
    float result = 0;
    int prim;
    prim = primByName(geo, vol, which);
    if (prim >= 0)
        result = volumesample(geo, prim, p);
    return result;
}

vector
volumesamplev(const string opname; const string prefix; const vector p; int which)
{
    string x, y, z;
    x = sprintf("%s.x", prefix);
    y = sprintf("%s.y", prefix);
    z = sprintf("%s.z", prefix);
    vector result = 0;
    int prim;
    prim = primByName(opname, x, which);
    if (prim >= 0)
        result.x = volumesample(opname, prim, p);
    prim = primByName(opname, y, which);
    if (prim >= 0)
        result.y = volumesample(opname, prim, p);
    prim = primByName(opname, z, which);
    if (prim >= 0)
        result.z = volumesample(opname, prim, p);
    return result;    
}

// Version of volumesamplev that allows selecting which vector volume to index.
vector
volumesamplev(const int geo; const string prefix; const vector p; int which)
{
    return volumesamplev(INPUT(geo), prefix, p, which);
}

float
amplitudeNoise(const string opname; const int prim;
               const float t; const vector p; const int plane)
{
    int noiseenabled = prim(opname, "noiseenabled", prim);
    if (!noiseenabled)
        return 1;

    // Fetch noise parameters from mask primitive.
    vector2 size = prim(opname, "noisesize", prim);
    int turb = prim(opname, "noiseturb", prim);
    float rough = prim(opname, "noiserough", prim);
    float dir = radians(prim(opname, "noisedir", prim));
    float pulsetime = prim(opname, "noisepulse", prim);
    float speed = prim(opname, "noisespeed", prim);
    vector offset = prim(opname, "noiseoffset", prim);
    int blend = prim(opname, "noiseblend", prim);
    vector2 input = prim(opname, "noiseinput", prim);
    vector2 output = prim(opname, "noiseoutput", prim);

    float scale = 1;

    // Vector in the direction of noise movement.
    vector noisevec = speed * threed(set(cos(dir), sin(dir)), plane);
    setplanar(noisevec, 1, plane);

    // Size and pulse time frequency.
    vector freq = threed(size, plane);
    setplanar(freq, pulsetime, plane);
    freq = 1 / freq;

    // Move and rotate in noise direction.
    matrix r = 1;
    vector up = 0;
    setplanar(up, 1, plane);

    translate(r, -noisevec * t);
    rotate(r, dir, up);
    scale(r, freq);
    translate(r, -offset*freq);

    // Turbulent sparse convolution noise.
    float n = 0;
    vector pp = p * r;
    for (int i = 0; i < turb; i++, pp *= 2, scale *= rough)
        n += scale * snoise(pp);

    // Map to ranges and use complement if blending.
    n = fit(n, input.x, input.y, output.x, output.y);
    if (blend)
        n = 1 - n;
    return n;
}

float maskAmplitude(const string maskname; const int layer;
                    const float t; const vector p; const int plane)
{
    int prim = primByName(maskname, "mask", layer);
    if (prim < 0)
        return 1;

    // Sample the mask in world-space.  Note that OceanSpectrum
    // SOP sets the border type to Streak, so 2-d masks should work
    // below the surface.
    float amp = volumesample(maskname, prim, p);
    amp = clamp(amp, 0, 1);
    // Suppression or contribution?
    int masktype = prim(maskname, "masktype", prim);
    if (masktype == 0)
        amp = 1 - amp;
    // Project input point back to the zero plane for noise,
    // so we get a consistent value no matter the depth.
    vector planarP = p;
    setplanar(planarP, 0, plane);
    amp *= amplitudeNoise(maskname, prim, t, planarP, plane);
    return amp;
}

float maskAmplitude(const int geo; const int layer;
                    const float t; const vector p; const int plane)
{
    return maskAmplitude(INPUT(geo), layer, t, p, plane);
}

int
hasPointMask(const string opname; const int layer)
{
    int prim = primByName(opname, "pointmask", layer);
    if (prim < 0)
        return 0;
    return len(primpoints(opname, prim)) > 0;
}

int
hasPointMask(const int geo; const int layer)
{
    return hasPointMask(INPUT(geo), layer);
}

void
antialiasDownsample(const vector2 res; const float gridsize; const vector p;
                    const float aablur; const int downsample;
                    int dsaa; float dsmix)
{
    // Calculate any downsampling required for anti-aliasing using the shading position.
    dsaa = 0;
    dsmix = 0;
    if (aablur > 0)
    {
        float minl = minwavelength(res, gridsize, downsample);
        float shadingarea = area(p);
        // Assume zero for very small values (e.g. from CVEX), 
        // so we get zero for dsmix and avoid a downsampled ocean_sample call.
        if (shadingarea < 1e-6)
            shadingarea = 0;
        else
            shadingarea = sqrt(shadingarea);
        float ratio = aablur * shadingarea / minl;
        dsaa = int(floor(ratio));
        // Maximum amount we can downsample.
        int maxds = int(floor(log(max(res)) / log(2))) - (downsample + 2);
        if (dsaa >= maxds)
        {
            dsaa = maxds;
            dsmix = 0;
        }
        else
        {
            dsmix = smooth(dsaa, dsaa + 1, ratio);
        }
    }
}

// Sample the actual ocean values.
void
sampleOceanValues(const string opname; const int isspectra;
                  const float maxdepth, time;
                  const vector p; const float ampscale, dsmix;
                  const int downsample, layer, plane;
                  vector displacement, v, dx)
{
    if (isspectra == 0)
    {
        // Already have evaluated rest volumes, just sample.
        // Clamp along planar axis if required when sampling in 3d.
        vector pclamp = p;
        if (maxdepth < 0)
            setplanar(pclamp, clamp(planar(pclamp, plane), maxdepth, 0), plane);
        displacement += ampscale * volumesamplev(opname, "restdisplace", pclamp, layer);
        v += ampscale * volumesamplev(opname, "restvel", pclamp, layer);
        dx += ampscale * volumesamplev(opname, "restderivatives", pclamp, layer);
        return;        
    }

    // Need to call ocean_sample.
    int phase = primByName(opname, "phase", layer);
    int freq = primByName(opname, "frequency", layer);
    int amp = primByName(opname, "amplitude", layer);
    float hscale = prim(opname, "chop", amp);
    float curtime = prim(opname, "timeoffset", amp);
    float timescale = 1;
    // Timescale attribute was recently introduced, check if exists first.
    if (hasprimattrib(opname, "timescale"))
        timescale = prim(opname, "timescale", amp);
    curtime += time * timescale;
    int ds = downsample;

    // Blend in first level's values.
    // We also need to scale v by any timescale.
    float scale = (1 - dsmix) * ampscale;
    displacement += scale * ocean_sample(opname, phase, freq, amp, hscale, curtime, 0, ds, p);
    v += timescale * scale * ocean_sample(opname, phase, freq, amp, hscale, curtime, 1, ds, p);
    dx += scale * ocean_sample(opname, phase, freq, amp, hscale, curtime, 2, ds, p);
    
    // Possibly blend in values downsampled by 1 more level.
    if (dsmix > 0)
    {
        scale = dsmix * ampscale;
        ds += 1;
        displacement += scale * ocean_sample(opname, phase, freq, amp, hscale, curtime, 0, ds, p);
        v += timescale * scale * ocean_sample(opname, phase, freq, amp, hscale, curtime, 1, ds, p);
        dx += scale * ocean_sample(opname, phase, freq, amp, hscale, curtime, 2, ds, p);
    }
}

// Derivative of smooth function, assuming min of 0.
float
dsmooth(const float p, r, z, rolloff)
{
    if (r >= p || r <= 0)
        return 0;
    float ds = 6 * (p - r) * z /  (p * p * p);
    if (rolloff > 1)
    {
        float h = r / p;
        float hs = h*h*(3 - 2*h);
        ds *= rolloff * pow(hs, rolloff - 1);
    }
    return ds;
}

// Sample a waveform of length "wavelength" attribute, with the wavefront
// oriented along the x-axis, spanning the z-axis.
void sampleWaveform(const string maskname; const int pt;
                    const vector pos; const matrix3 Rinv;
                    const float amporig, S, pscale, r, rolloff, depthfalloff;
                    const int plane;
                    vector displacement, v, dx)
{
    float l = point(maskname, "wavelength", pt);
    float halfwave = l / 2;
    float xa = abs(pos.x);
    // Nothing to do outside wavelength.
    // This will always exit on zero (or non-existent) wavelength.
    if (xa >= halfwave)
        return;
    float amp = amporig;
    // k is wavenumber.
    float k = abs(2 * M_PI / l);
    // Negative falloff means exponential by frequency
    if (depthfalloff < 0) 
        amp *= exp(-depthfalloff * k * min(0, planar(pos, plane)));
    else if (depthfalloff > 0) 
        amp *= exp(depthfalloff * min(0, planar(pos, plane)));
    float cosw = cos(k * pos.x);
    float sinw = sin(k * pos.x);
    // From "Wave Particles", Yuksel et al, SIGGRAPH 2007
    float D = 0.5 * amp * (cosw + 1);
    // Y displacement is waveform * smoothing function.
    displacement.y += D * S;

    // Gerstner wave for X displacement
    float X = -amp * sinw;
    float dXdx = (-k * amp) * cosw;
    // Smooth horizontal falloff for better wave shape and zero
    // first derivative at limits.
    float Br = 10;
    float B = 1 - smooth(0, halfwave, xa, Br);
    float dBdx = -dsmooth(halfwave, xa, pos.x, Br);
    
    // X displacement is chop * horizontal waveform * smoothing functions.
    float chop = point(maskname, "wavechop", pt);
    displacement.x += chop * X * S * B;

    // Spatial derivatives of S in x.
    vector offset = point(maskname, "offset", pt);
    // We calculated r in world space, so need to take offset into account
    // when computing derivatives at the current position (also for dSdz below).
    float dSdx = -dsmooth(pscale, r, pos.x - offset.x, rolloff);

    // Spatial derivative of x += chop * X * S * B via product rule:
    dx.x += chop * (X * dSdx * B  + S * dXdx * B + X * S * dBdx);

    // Add in wave velocity if non-zero point velocity.
    vector pv = point(maskname, "v", pt);
    if (length2(pv) > 0)
    {
        // Rotate velocity to waveform space.
        pv = pv * Rinv;
        // Time derivatives of x and z are negative local v.
        float dxdt = -pv.x;
        float dzdt = -pv.z;
        // Time derivative of D in x (invariant in z).
        float dD_xdx = (-M_PI * amp / l) * sinw;
        float dD_xdt = dD_xdx * dxdt;

        // Spatial derivative of S in z.
        float dSdz = -dsmooth(pscale, r, pos.z - offset.z, rolloff);
        // Time derivative of S in x and z directions.
        float dS_xdt = dSdx * dxdt;
        float dS_zdt = dSdz * dzdt;
        // Time derivative of y += D * S via product rule:
        v.y += D * (dS_xdt + dS_zdt) + S * dD_xdt;

        // Time derivative of L and B in x (invariant in z).
        float dX_xdt = dXdx * dxdt;
        float dB_xdt = dBdx * dxdt;
        // Time derivative of x += chop * L * S * B via product rule.
        v.x += chop * (X * (dS_xdt + dS_zdt) * B + S * dX_xdt * B + X * S * dB_xdt);
    }
}

int
accumulatePointInstances(const string opname; const string maskname;
                         const int isspectra; const float maxdepth, time; 
                         const vector p; const float maskamp, dsmix, depthfalloff;
                         const int downsample, layer, plane;
                         vector displacement, v; matrix3 J)
{
    int ptmask = primByName(maskname, "pointmask", layer);
    if (ptmask < 0)
        return 0;
    int nvtx = primvertexcount(maskname, ptmask);
    if (nvtx == 0)
        return 0;

    // Variable radius point lookup.
    // First project lookup point to best-fit plane, giving a
    // cylindrical lookup for proper velocities at depth.
    vector pplanar = p;
    vector planeN = prim(maskname, "planeN", ptmask);
    if (length2(planeN) > 0)
    {
        // Project lookup point to best-fit plane.
        vector planeorig = prim(maskname, "planeorig", ptmask);
        pplanar -= dot(pplanar - planeorig, planeN) * planeN;
    }

    // Create a point group pattern specifying pointmask.
    // NOTE:  This assumes that the vertex order of the points is sorted
    // ascending, which *should* always be the case if the particle system primitive
    // was created by OceanSpectrum SOP. This will break if the user shuffles
    // the order of the wave instance points after OceanSpectrum, so...
    // they shouldn't do that.
    int minpt = vertexpoint(maskname, primvertex(maskname, ptmask, 0));
    int maxpt = vertexpoint(maskname, primvertex(maskname, ptmask, nvtx-1));
    string ptgrp = sprintf("%d-%d", minpt, maxpt);
    int pts[] = pcfind_radius(maskname, ptgrp, "P", "pscale", 1, pplanar, 0, 100000);

    // Check if any points have waveform attribute, and if we need to sample spectral values.
    int hasspectral = 1;
    int haswaveform = haspointattrib(maskname, "wavelength");
    if (haswaveform && plane == 2)
    {
        int ampprim = -1;
        // Find the representative volume for the spectrum that determines if null amplitudes.
        if (isspectra)
        {
            ampprim = primByName(maskname, "amplitude", layer);
        }
        else
        {
            ampprim = primByName(maskname, "restdisplace.y", layer);
            if (ampprim < 0)
                ampprim = primByName(maskname, "restvel.y", layer);
            if (ampprim < 0)
                ampprim = primByName(maskname, "restderivatives.y", layer);
        }
        // We're specifically trying to avoid sampling from the y-up layer created by the OceanWaves SOP,
        // which has null spectral volumes of res 128.  Since the min/max volume intrinsics are not
        // cached, we only check when it's cheap to do so.
        if (ampprim >= 0 && max(volumeres(maskname, ampprim)) <= 128)
            hasspectral = (abs(primintrinsic(maskname, "volumemaxvalue", ampprim)) > 1e-6 ||
                           abs(primintrinsic(maskname, "volumeminvalue", ampprim)) > 1e-6);
    }

    float rolloff = prim(maskname, "rolloff", ptmask);
    foreach(int pt; pts)
    {
        float r = distance(pplanar, point(maskname, "P", pt));
        float pscale = point(maskname, "pscale", pt);
        if (r > pscale)
            continue;
        // Get point instance transform inverse.
        matrix xforminv = point(maskname, "xforminv", pt);
        
        // Transform back to ocean space.
        vector pos = p * xforminv;
        // Rotation matrix for transforming sampled vectors.
        matrix3 Rinv = matrix3(xforminv);
        matrix3 R = transpose(Rinv);

        // Scale the mask amplitude by per-point contributions.
        float amp = maskamp * point(maskname, "amplitude", pt);
        // Smoothing function for wave amplitude falloff.
        float S = 1 - smooth(0, pscale, r, rolloff);
        
        // Accumulate ocean-space displacement, velocity, and derivatives.
        vector ldisp = 0, lv = 0, ldx = 0;
        if (abs(amp * S) > 0)
        {
            // Sample from the ocean spectrum in ocean space.
            if (hasspectral)
                sampleOceanValues(opname, isspectra, maxdepth, time, pos, amp * S, dsmix,
                                  downsample, layer, plane, ldisp, lv, ldx);

            // Sample from a single waveform or length "wavelength".
            if (haswaveform)
                sampleWaveform(maskname, pt, pos, Rinv, amp, S, pscale, r,
                               rolloff, depthfalloff, plane,
                               ldisp, lv, ldx);
        }

        // Rotate and accumulate world-space displacements and velocities.
        displacement += ldisp * R;
        v += lv * R;
        // Calc Jacobian of rotated derivatives from inverse rotated lookup position.
        // J(D(pos * Rinv) * R) =
        // J(Rinv) * J(D(pos)) * J(R) = 
        // Rinv * J(D(pos)) * R
        matrix3 lJ = jacobian3d(ldx, plane);
        J += Rinv * lJ * R;
    }
    // Assume success if any points exist, even if none contribute to displacement.
    return 1;
}

void
oceanSampleLayers(const string opname; const string maskname; 
            const float time; const vector p;
            const float aablur; const float depthfalloff; const int downsample;
            vector displacement; vector v; matrix3 J)
{
    // Count layers via either spectrum or evaluated volume names.
    string volumes[] = { "amplitude", "restdisplace.x", "restvel.x", "restderivatives.x" };
    int nlayers = 0;
    string vol;
    // Assume we're getting evaluated volumes.
    int isspectra = 0;
    for(int i = 0; i < len(volumes); i++)
    {
        vol = volumes[i];
        nlayers = nprimsByName(opname, vol);
        if (nlayers)
        {
            // Spectrum or evaluated rest input?
            isspectra = (vol == "amplitude");
            break;
        }
    }

    // Check for valid separate mask input, which should have
    // a "mask" and "pointmask" primitive per layer.
    int hasmask = (nprimsByName(maskname, "mask") == nlayers) &&
                  (nprimsByName(maskname, "pointmask") == nlayers);

    displacement = 0;
    v = 0;
    J = 0;

    // Check if nothing to do.
    if (!nlayers)
        return;

    for(int i=0; i < nlayers; i++)
    {
        int prim = primByName(opname, vol, i);
        // Resolution of current layer.
        vector res3 = volumeres(opname, prim);
        int plane = plane2d(res3);
        vector2 res = twod(res3, plane);

        // Grid size of current layer.
        string primgrp = sprintf("%d", prim);
        vector bbmin, bbmax;
        getbbox(opname, primgrp, bbmin, bbmax);
        vector size = bbmax - bbmin;
        float gridsize = max(size);

        // Calc maxdepth for clamping before sampling from evaluated volumes.
        float maxdepth = 0;
        float resp = planar(res3, plane);
        // Max sample depth is 1/2 voxel from bottom.
        if (!isspectra && resp > 1)
            maxdepth = planar(bbmin, plane) + planar(size, plane) / (2 * resp);

        // Mask scaling for this layer.
        float amp = maskAmplitude(hasmask ? maskname : opname, i, time, p, plane);

        // Calculate downsample values for anti-aliasing at this shading point.
        int dsaa;
        float dsmix;
        antialiasDownsample(res, gridsize, p, aablur, downsample, dsaa, dsmix);

        // If point instances contribute, we're done with this layer.
        if(accumulatePointInstances(opname, hasmask ? maskname : opname,
                                    isspectra, maxdepth, time, p, amp,
                                    dsmix, depthfalloff, downsample + dsaa,
                                    i, plane,
                                    displacement, v, J))
            continue;

        // Otherwise sample from this position as usual.
        vector dx = 0;
        sampleOceanValues(opname, isspectra, maxdepth, time, p, amp,
                          dsmix, downsample + dsaa, i, plane,
                          displacement, v, dx);
        J += jacobian3d(dx, plane);
    }
}

// Versions with no anti-aliasing, no downsampling.
void
oceanSampleLayers(const string opname; const float time; const vector p;
                  const float depthfalloff;
                  vector displacement; vector v; matrix3 J)
{
    oceanSampleLayers(opname, "", time, p, 0, depthfalloff, 0, displacement, v, J);
}

void
oceanSampleLayers(const int geo; const float time; const vector p;
                  const float depthfalloff;
                  vector displacement; vector v; matrix3 J)
{
    oceanSampleLayers(INPUT(geo), "", time, p, 0, depthfalloff, 0, displacement, v, J);
}

// Version that compute cusp values directly, called by OceanSampleLayers VOP.
void
oceanSampleLayers(const string opname; const string maskname; const float time; const vector p;
                  const float aablur; const int falloffmode; const float falloffscale; const int downsample;
                  vector displacement; vector v; float cusp; vector cuspdir)
{
    matrix3 J;
    float falloff = 0;
    if (falloffmode == 1)
        falloff = falloffscale;
    else if (falloffmode == 2)
        falloff = -falloffscale;
    oceanSampleLayers(opname, maskname, time, p, aablur, falloff, downsample, displacement, v, J);
    computeCusp(J, cusp, cuspdir, oceanPlane2d(opname));
}

#endif
