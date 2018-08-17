/****************************************************************************
   Description: Functions to generate samples with various 2D/3D
                distributions. Used by pyro to distribute samples
                inside the differential space of a volume position.

*****************************************************************************/

#ifndef pyro_sampler_h__GUARD
#define pyro_sampler_h__GUARD

#include <pyro_math.h>

// Uniforms for the gaussian inverse cdf approximation
#define  G_A1  -3.969683028665376e+01
#define  G_A2   2.209460984245205e+02
#define  G_A3  -2.759285104469687e+02
#define  G_A4   1.383577518672690e+02
#define  G_A5  -3.066479806614716e+01
#define  G_A6   2.506628277459239e+00

#define  G_B1  -5.447609879822406e+01
#define  G_B2   1.615858368580409e+02
#define  G_B3  -1.556989798598866e+02
#define  G_B4   6.680131188771972e+01
#define  G_B5  -1.328068155288572e+01

#define  G_C1  -7.784894002430293e-03
#define  G_C2  -3.223964580411365e-01
#define  G_C3  -2.400758277161838e+00
#define  G_C4  -2.549732539343734e+00
#define  G_C5   4.374664141464968e+00
#define  G_C6   2.938163982698783e+00

#define  G_D1   7.784695709041462e-03
#define  G_D2   3.224671290700398e-01
#define  G_D3   2.445134137142996e+00
#define  G_D4   3.754408661907416e+00

#define  G_P_LOW   0.02425
#define  G_P_HIGH  0.97575    // P_HIGH = 1 - p_LOW
#define  G_W       0.999967   // indexing limit to the gaussian icdf



//-------------------------------------------------------------------------------
// Filter Distributions
// These functions represent random variates in N dimensions (N in {1,2,3}) of the
// filter functions supported by pyro2. They take N uniformly distributed random
// numbers and return a position vector (of length N) with the chosen probability 
// distribution. They also export the pdf associated with the sample.
// These filters are elliptical (anisotropic) -- they can take different
// widths in each dimension.
// The naming convension for these functions is:
//   <base_name><dim>_<filter_type>, where:
//   base_name   = "fsample" (for "filter sample")
//   dim         = one of: {1d, 2d, 3d}
//   filter_type = one of: {box, bartlett, gauss}
// Assumptions:  The input random numbers are in [0,1) or [0,1]
// Ouput Ranges: fscale*[-1,1] or fscale*[-2,2] depending on the filter type
//-------------------------------------------------------------------------------

// Concentric sampling of disk
// this version returns the components (r,sintheta,costheta) which are
// used by some of the other sampling functions
void pdist_diskcomps(float u1, u2; export float r,th,cth,sth) 
{
   r = th = cth = 0;
   sth = 1;

   // Map uniform random numbers to [-1,1]^2
   float sx = 2 * u1 - 1;
   float sy = 2 * u2 - 1;
   
   // Handle degeneracy at the origin
   if(sx==0.0 && sy==0.0) return;
   
   if(sx >= -sy) {
       if(sx > sy) {
           // Handle first region of disk
           r = sx;
           if(sy > 0) th = sy/r;
           else       th = 8 + sy/r;
       } else {
           // Handle second region of disk
           r = sy;
           th = 2 - sx/r;
       }
   } else {
       if(sx <= sy) {
           // Handle third region of disk
           r = -sx;
           th = 4 - sy/r;
       } else {
           // Handle fourth region of disk
           r = -sy;
           th = 6 + sx/r;
       }
   }
   th *= C_PI_4;
   cth = cos(th);
   sth = sin(th);
}

// this version simply returns the final position
vector sdist_disk(float u1, u2) {
   float r,th,cth,sth;
   pdist_diskcomps(u1,u2, r,th,cth,sth);
   return r * set(cth, sth, 0);
}
vector sdist_disk(vector uv) {
   return sdist_disk(uv[0],uv[1]);
}


//------------------------------------------------------------------------------
// Sphere::Uniform
//------------------------------------------------------------------------------
void pdist_spherecomps(float u1,u2; float up,r,ph) {
   up = 1 - 2*u1;
   r  = sqrt(max(0, 1 - up*up));
   ph = C_PI2*u2;
}
// Orientation: +x
vector pdist(Sphere dom; Uniform pdf; XAxis orient; float u1,u2) {
   float x,r,ph; pdist_spherecomps(u2,u1, x,r,ph);
   return set(-x, r*sin(ph), -r*cos(ph));
}
vector pdist(Sphere dom; Uniform pdf; XAxis orient; vector u) {
   return pow(u.z,C_1_3) * pdist(dom,pdf,orient, u.x,u.y);
}
// Orientation: +y
vector pdist(Sphere dom; Uniform pdf; YAxis orient; float u1,u2) {
   float y,r,ph; pdist_spherecomps(u2,u1, y,r,ph);
   return set(-r*cos(ph), -y, r*sin(ph));
}
vector pdist(Sphere dom; Uniform pdf; YAxis orient; vector u) {
   return pow(u.z,C_1_3) * pdist(dom,pdf,orient, u.x,u.y);
}
// Orientation: +z
vector pdist(Sphere dom; Uniform pdf; ZAxis orient; float u2,u1) {
   float z,r,ph; pdist_spherecomps(u1,u2, z,r,ph);
   return -set(r*cos(ph), r*sin(ph), z);
}
vector pdist(Sphere dom; Uniform pdf; ZAxis orient; vector u) {
   return pow(u.z,C_1_3) * pdist(dom,pdf,orient, u.x,u.y);
}
// Default Orientation: +z
vector pdist(Sphere dom; Uniform pdf; float u2,u1) {
   return pdist(dom,pdf,ZAxis(), u1,u2);
}
vector pdist(Sphere dom; Uniform pdf; vector u) {
   return pdist(dom,pdf,ZAxis(), u);
}
// Default PDF: Uniform
vector pdist(Sphere dom; XAxis orient; float u2,u1) {
   return pdist(dom,Uniform(),orient, u1,u2);
}
vector pdist(Sphere dom; XAxis orient; vector u) {
   return pdist(dom,Uniform(),orient, u);
}
vector pdist(Sphere dom; YAxis orient; float u2,u1) {
   return pdist(dom,Uniform(),orient, u1,u2);
}
vector pdist(Sphere dom; YAxis orient; vector u) {
   return pdist(dom,Uniform(),orient, u);
}
vector pdist(Sphere dom; ZAxis orient; float u2,u1) {
   return pdist(dom,Uniform(),orient, u1,u2);
}
vector pdist(Sphere dom; ZAxis orient; vector u) {
   return pdist(dom,Uniform(),orient, u);
}
vector pdist(Sphere dom; float u2,u1) {
   return pdist(dom,Uniform(),ZAxis(), u1,u2);
}
vector pdist(Sphere dom; vector u) {
   return pdist(dom,Uniform(),ZAxis(), u);
}
// Default Domain: Sphere
vector pdist(XAxis orient; float u2,u1) {
   return pdist(Sphere(),Uniform(),orient, u1,u2);
}
vector pdist(XAxis orient; vector u) {
   return pdist(Sphere(),Uniform(),orient, u);
}
vector pdist(YAxis orient; float u2,u1) {
   return pdist(Sphere(),Uniform(),orient, u1,u2);
}
vector pdist(YAxis orient; vector u) {
   return pdist(Sphere(),Uniform(),orient, u);
}
vector pdist(ZAxis orient; float u2,u1) {
   return pdist(Sphere(),Uniform(),orient, u1,u2);
}
vector pdist(ZAxis orient; vector u) {
   return pdist(Sphere(),Uniform(),orient, u);
}
vector pdist(float u2,u1) {
   return pdist(Sphere(),Uniform(),ZAxis(), u1,u2);
}
vector pdist(vector u) {
   return pdist(Sphere(),Uniform(),ZAxis(), u);
}


//------------------------------------------------------------------------------
// Hemisphere::Uniform
//------------------------------------------------------------------------------
void pdist_hemicomps(float u1,u2, up,r,ph) {
   up = u1;
   r  = sqrt(max(0, 1 - up*up));
   ph = C_PI2 * u2;
}
// Orientation: +x
vector pdist(Hemisphere dom; Uniform pdf; XAxis orient; float u1,u2) {
   float x,r,ph; pdist_hemicomps(u2,u1, x,r,ph);
   return set(x, r*sin(ph), -r*cos(ph));
}
vector pdist(Hemisphere dom; Uniform pdf; XAxis orient; vector u) {
   return pow(u.z,C_1_3) * pdist(dom,pdf,orient, u.x,u.y);
}
// Orientation: +y
vector pdist(Hemisphere dom; Uniform pdf; YAxis orient; float u1,u2) {
   float y,r,ph; pdist_hemicomps(u2,u1, y,r,ph);
   return set(-r*cos(ph), y, r*sin(ph));
}
vector pdist(Hemisphere dom; Uniform pdf; YAxis orient; vector u) {
   return pow(u.z,C_1_3) * pdist(dom,pdf,orient, u.x,u.y);
}
// Orientation: +z
vector pdist(Hemisphere dom; Uniform pdf; ZAxis orient; float u1,u2) {
   float z,r,ph; pdist_hemicomps(u2,u1, z,r,ph);
   return set(-r*cos(ph), -r*sin(ph), z);
}
vector pdist(Hemisphere dom; Uniform pdf; ZAxis orient; vector u) {
   return pow(u.z,C_1_3) * pdist(dom,pdf,orient, u.x,u.y);
}
// Default Orientation: +z
vector pdist(Hemisphere dom; Uniform pdf; float u2,u1) {
   return pdist(dom,pdf,ZAxis(), u1,u2);
}
vector pdist(Hemisphere dom; Uniform pdf; vector u) {
   return pdist(dom,pdf,ZAxis(), u);
}
// Default PDF: Uniform
vector pdist(Hemisphere dom; XAxis orient; float u2,u1) {
   return pdist(dom,Uniform(),orient, u1,u2);
}
vector pdist(Hemisphere dom; XAxis orient; vector u) {
   return pdist(dom,Uniform(),orient, u);
}
vector pdist(Hemisphere dom; YAxis orient; float u2,u1) {
   return pdist(dom,Uniform(),orient, u1,u2);
}
vector pdist(Hemisphere dom; YAxis orient; vector u) {
   return pdist(dom,Uniform(),orient, u);
}
vector pdist(Hemisphere dom; ZAxis orient; float u2,u1) {
   return pdist(dom,Uniform(),orient, u1,u2);
}
vector pdist(Hemisphere dom; ZAxis orient; vector u) {
   return pdist(dom,Uniform(),orient, u);
}
vector pdist(Hemisphere dom; float u2,u1) {
   return pdist(dom,Uniform(),ZAxis(), u1,u2);
}
vector pdist(Hemisphere dom; vector u) {
   return pdist(dom,Uniform(),ZAxis(), u);
}


// A fast algorithm (by Peter Acklam, http://home.online.no/~pjacklam/notes/invnorm)
// for computing the inverse cdf of a normal distribution with stdev=1 and mean=0.
// Argument "u" is meant to be a uniformly distributed random number in [0,1]

float icdf_gauss(float u, mean, stdev) {
   // Note: Acklam's algorithm was built around the canonical gaussian centered
   // at 0 with stdev=1, because the other parameterizations are just scalings
   // and offsets of this. It also (correctly) has an [-inf,+inf] ouput range. 
   // But, we'll be using it for drawing samples from the finite segment that
   // is the width of the filter, so we need to limit its range. This is 
   // typically done by windowing the PDF, forcing the tails to reach 0 at some
   // specific radius. Here, we'll limit the indexing to the icdf instead.
   // The following line squeezes the input index (in [0,1]) such that the
   // icdf has an output range of [-2,2] (which is the customary window size
   // for the standard gaussian). It is the only change to Acklam's algorithm.
   float v = fit(u,0,1,G_W,1-G_W);

   float x=0;
   float q,r;
   if((0 < v )  && (v < G_P_LOW)) {
      q = sqrt(-2*log(v));
      x = (((((G_C1*q+G_C2)*q+G_C3)*q+G_C4)*q+G_C5)*q+G_C6) / 
          ((((G_D1*q+G_D2)*q+G_D3)*q+G_D4)*q+1);
   } else {
      if((G_P_LOW <= v) && (v <= G_P_HIGH)) {
         q = v - 0.5;
         r = q*q;
         x = (((((G_A1*r+G_A2)*r+G_A3)*r+G_A4)*r+G_A5)*r+G_A6)*q /
             (((((G_B1*r+G_B2)*r+G_B3)*r+G_B4)*r+G_B5)*r+1);
      } else {
         if((G_P_HIGH < v) && (v < 1)) {
            q = sqrt(-2*log(1-v));
            x = -(((((G_C1*q+G_C2)*q+G_C3)*q+G_C4)*q+G_C5)*q+G_C6) / 
                 ((((G_D1*q+G_D2)*q+G_D3)*q+G_D4)*q+1);
         }
      }
   }
   return stdev * -x + mean;
}

float pdf_gauss(float x, mean, stdev) {
   float u  = x - mean;
   float s2 = stdev*stdev;
   return exp(-(u*u) / (2*s2)) / (C_SQRT2PI * stdev);
}

float icdf_bartlett(float u, mean) {
   float x = (u>0.5) ? 1 - sqrt(2 - 2*u) : C_SQRT2*sqrt(u) - 1;
   return 0.5 * x + mean;
}

float pdf_bartlett(float x, mean) {
   float  u = abs(x-mean);
   return 2 * max(0, 1 - 2*u);
}


// Mean and standard deviation for the Gaussian distributions here.
// They're set so that its natural range is [-1,1] to match the 
// uniform extents above
#define PDGAUSS_MEAN    0     
#define PDGAUSS_STDEV   0.25


//------------------------------------------------------------------------------
// Sphere::Gaussian
//------------------------------------------------------------------------------
// 1D
float pdist(float u; Sphere dom; Gaussian dist; XAxis o) {
   return icdf_gauss(u,PDGAUSS_MEAN,PDGAUSS_STDEV);
}
float pdist(Sphere dom; Gaussian dist; XAxis o; float u; export float pdf) {
   float x = pdist(u, dom,dist,o);
   pdf = pdf_gauss(x,PDGAUSS_MEAN,PDGAUSS_STDEV);
   return x;
}
// 2D
vector pdist(Sphere dom; Gaussian dist; XAxis o; float ux,uy) {
   float r,th,cth,sth;
   pdist_diskcomps(ux,uy, r,th,cth,sth);
   return icdf_gauss(r*.5+.5,PDGAUSS_MEAN,PDGAUSS_STDEV) * set(cth,sth,0);
}
vector pdist(Sphere dom; Gaussian dist; XAxis o; float ux,uy; export float pdf) {
   float r,th,cth,sth;
   pdist_diskcomps(ux,uy, r,th,cth,sth);
   r = icdf_gauss(r*.5+.5,PDGAUSS_MEAN,PDGAUSS_STDEV);
   pdf = pdf_gauss(r,PDGAUSS_MEAN,PDGAUSS_STDEV);
   return r*set(cth,sth,0);
}
// 3D
vector pdist(Sphere dom; Gaussian dist; XAxis o; vector u) {
   vector dir = pdist(dom,dist,o, u.x,u.y);
   return dir * icdf_gauss(pow(u.z,C_1_3)*.5+.5,PDGAUSS_MEAN,PDGAUSS_STDEV);
}
vector pdist(Sphere dom; Gaussian dist; XAxis o; vector u; export float pdf) {
   vector dir = pdist(dom,dist,o, u.x,u.y);
   float  r   = icdf_gauss(pow(u.z,C_1_3)*.5+.5,PDGAUSS_MEAN,PDGAUSS_STDEV);
   pdf = pdf_gauss(r,PDGAUSS_MEAN,PDGAUSS_STDEV);
   return r * dir;
}


// Bartlett (triangle) Filter ( width: 1, output range: [-0.5,0.5] )
//--------------------------------------------------------------------------
// 1D
float  sdist_bartlett(float u) {
   return icdf_bartlett(u,0);
}
float  sdist_bartlett(float u; export float pdf) {
   float x = sdist_bartlett(u);
   pdf = pdf_bartlett(x,0);
   return x;
}
// 2D
vector sdist_bartlett(float ux,uy) {
   float r,th,cth,sth; pdist_diskcomps(ux,uy, r,th,cth,sth);
   r = icdf_bartlett(r,0);
   return r * set(cth,sth,0);
}
vector sdist_bartlett(float ux,uy; export float pdf) {
   float r,th,cth,sth; pdist_diskcomps(ux,uy, r,th,cth,sth);
   r = icdf_bartlett(r,0);
   pdf = pdf_bartlett(r,0);
   return r * set(cth,sth,0);
}

#endif // End pyro_sampler_h
