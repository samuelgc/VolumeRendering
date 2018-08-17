/****************************************************************************
   Description:  Basic utils for pyro.

*****************************************************************************/

#ifndef pyro_utils_h__GUARD
#define pyro_utils_h__GUARD

#include <pyro_math.h>

#define CAT2_IMPL(a,b)		a##b
#define CAT2(a,b)     		CAT2_IMPL(a,b)
#define CAT3_IMPL(a,b,c)	a##b##c
#define CAT3(a,b,c)   		CAT3_IMPL(a,b,c)
#define CAT4_IMPL(a,b,c,d) 	a##b##c##d
#define CAT4(a,b,c,d) 		CAT4_IMPL(a,b,c,d)


#define STRINGIFY_IMPL(a)  #a
#define STRINGIFY(a)  STRINGIFY_IMPL(a)

//-------------------------------------------------------------------------------
// The VOP Context (one of SOP,POP,COP,CHOP,SHADING)
//-------------------------------------------------------------------------------
#if defined(VOP_CTXT)
   #undef VOP_CTXT
#endif

#if defined(VOP_SHADING)

   #define VOP_CTXT  SHADING
   #if defined(VOP_SURFACE)
      #define VOP_SHCTXT SURFACE
   #elif defined(VOP_LIGHT)
      #define VOP_SHCTXT LIGHT
   #elif defined(VOP_SHADOW)
      #define VOP_SHCTXT SHADOW
   #elif defined(VOP_FOG)
      #define VOP_SHCTXT FOG
   #elif defined(VOP_DISPLACE)
      #define VOP_SHCTXT DISPLACEMENT
   #endif

#elif defined(VOP_SOP)
   #define VOP_CTXT  SOP

#elif defined(VOP_POP)
   #define VOP_CTXT  POP

#elif defined(VOP_COP) || defined(VOP_COP2)
   #define VOP_CTXT  COP

#elif defined(VOP_CHOP)
   #define VOP_CTXT  CHOP

#endif


//-------------------------------------------------------------------------------
// The type identifiers used for function names, just to enforce some consistency
//-------------------------------------------------------------------------------
#define TYPEIDfloat    f
#define TYPEIDvector   v
#define TYPEIDvector4  v4
#define TYPEIDmatrix3  m3
#define TYPEIDmatrix   m



//-------------------------------------------------------------------------------
// Wave Forms.
// Range: [0,1]
// Default: Sine, period 1, unfiltered
//-------------------------------------------------------------------------------

// Triangle, Unfiltered
float wave(Triangle type; float x) {
   return 1 - 2*abs(0.5-modulo(x,1));
}
float wave(Triangle type; float x, period) {
   float per = exceptzero(period,1e-5);
   return 1 - 2*abs(0.5-modulo(x,per)/per);
}
vector wave(Triangle type; vector x) {
   return VONE - 2*abs((vector)0.5-modulo(x,1));
}
vector wave(Triangle type; vector x, period) {
   vector per = exceptzero(period,1e-5);
   return VONE - 2*abs((vector)0.5-modulo(x,per)/per);
}
vector wave(Triangle type; vector x; float period) {
   float per = exceptzero(period,1e-5);
   return VONE - 2*abs((vector)0.5-modulo(x,per)/per);
}
vector4 wave(Triangle type; vector4 x) {
   return (vector4)1 - 2*abs((vector4)0.5-modulo(x,1));
}
vector4 wave(Triangle type; vector4 x, period) {
   vector4 per = exceptzero(period,1e-5);
   return (vector4)1 - 2*abs((vector4)0.5-modulo(x,per)/per);
}
vector4 wave(Triangle type; vector4 x; float period) {
   float per = exceptzero(period,1e-5);
   return (vector4)1 - 2*abs((vector4)0.5-modulo(x,per)/per);
}

// Cosine, Unfiltered
float wave(Cosine type; float x) {
   float y = float(wave(Triangle(),x));
   return cos(y*C_PI)*.5+.5;
}
float wave(Cosine type; float x, per) {
   float y = float(wave(Triangle(),x,per));
   return cos(y*C_PI)*.5+.5;
}
vector wave(Cosine type; vector x) {
   vector y = vector(wave(Triangle(),x));
   return cos(y*C_PI)*.5+.5;
}
vector wave(Cosine type; vector x, per) {
   vector y = vector(wave(Triangle(),x,per));
   return cos(y*C_PI)*.5+.5;
}
vector wave(Cosine type; vector x; float per) {
   vector y = vector(wave(Triangle(),x,per));
   return cos(y*C_PI)*.5+.5;
}
vector4 wave(Cosine type; vector4 x) {
   vector4 y = vector4(wave(Triangle(),x));
   return cos(y*C_PI)*.5+.5;
}
vector4 wave(Cosine type; vector4 x, per) {
   vector4 y = vector4(wave(Triangle(),x,per));
   return cos(y*C_PI)*.5+.5;
}
vector4 wave(Cosine type; vector4 x; float per) {
   vector4 y = vector4(wave(Triangle(),x,per));
   return cos(y*C_PI)*.5+.5;
}

// Sine, Unfiltered
float wave(Sine type; float x) {
   return 1 - float(wave(Cosine(),x));
}
float wave(Sine type; float x, per) {
   return 1 - float(wave(Cosine(),x,per));
}
vector wave(Sine type; vector x) {
   return VONE - vector(wave(Cosine(),x));
}
vector wave(Sine type; vector x, per) {
   return VONE - vector(wave(Cosine(),x,per));
}
vector wave(Sine type; vector x; float per) {
   return VONE - vector(wave(Cosine(),x,per));
}
vector4 wave(Sine type; vector4 x) {
   return (vector4)1 - vector4(wave(Cosine(),x));
}
vector4 wave(Sine type; vector4 x, per) {
   return (vector4)1 - vector4(wave(Cosine(),x,per));
}
vector4 wave(Sine type; vector4 x; float per) {
   return (vector4)1 - vector4(wave(Cosine(),x,per));
}

// Cubic (Smooth), Unfiltered
float wave(Cubic type; float x) {
   float y = float(wave(Triangle(),x));
   return (3 - 2*y)*y*y;
}
float wave(Cubic type; float x, per) {
   float y = float(wave(Triangle(),x,per));
   return (3 - 2*y)*y*y;
}
vector wave(Cubic type; vector x) {
   vector y = vector(wave(Triangle(),x));
   return ((vector)3 - 2*y)*y*y;
}
vector wave(Cubic type; vector x, per) {
   vector y = vector(wave(Triangle(),x,per));
   return ((vector)3 - 2*y)*y*y;
}
vector wave(Cubic type; vector x; float per) {
   vector y = vector(wave(Triangle(),x,per));
   return ((vector)3 - 2*y)*y*y;
}
vector4 wave(Cubic type; vector4 x) {
   vector4 y = vector4(wave(Triangle(),x));
   return ((vector4)3 - 2*y)*y*y;
}
vector4 wave(Cubic type; vector4 x, per) {
   vector4 y = vector4(wave(Triangle(),x,per));
   return ((vector4)3 - 2*y)*y*y;
}
vector4 wave(Cubic type; vector4 x; float per) {
   vector4 y = vector4(wave(Triangle(),x,per));
   return ((vector4)3 - 2*y)*y*y;
}

// Sawtooth, Unfiltered
float wave(Sawtooth type; float x) {
   return modulo(x,1);
}
float wave(Sawtooth type; float x, period) {
   float per = exceptzero(period,1e-5);
   return modulo(x,per)/per;
}
vector wave(Sawtooth type; vector x) {
   return modulo(x,1);
}
vector wave(Sawtooth type; vector x, period) {
   vector per = exceptzero(period,1e-5);
   return modulo(x,per)/per;
}
vector wave(Sawtooth type; vector x; float period) {
   float per = exceptzero(period,1e-5);
   return modulo(x,per)/per;
}
vector4 wave(Sawtooth type; vector4 x) {
   return modulo(x,1);
}
vector4 wave(Sawtooth type; vector4 x, period) {
   vector4 per = exceptzero(period,1e-5);
   return modulo(x,per)/per;
}
vector4 wave(Sawtooth type; vector4 x; float period) {
   float per = exceptzero(period,1e-5);
   return modulo(x,per)/per;
}


// Sine is the default wave type
float wave(float x) {
   return float(wave(Sine(),x));
}
float wave(float x, per) {
   return float(wave(Sine(),x,per));
}
vector wave(vector x) {
   return vector(wave(Sine(),x));
}
vector wave(vector x, per) {
   return vector(wave(Sine(),x,per));
}
vector wave(vector x; float per) {
   return vector(wave(Sine(),x,per));
}
vector4 wave(vector4 x) {
   return vector4(wave(Sine(),x));
}
vector4 wave(vector4 x, per) {
   return vector4(wave(Sine(),x,per));
}
vector4 wave(vector4 x; float per) {
   return vector4(wave(Sine(),x,per));
}



//-------------------------------------------------------------------------------
// Arrays
//-------------------------------------------------------------------------------

// Integrate an array.
// Assumes array entries are samples taken on even steps of length 'xstep'
// and does trapezoidal integration over the entire array.
float ainteg(float v[], xstep) {
   int   len = len(v), i;
   float fa  = v[0];
   if(len<2) return fa;
   float out = 0, fb;

   for(i=1; i<len; i++) {
      fb = v[i];
      out += fa + fb;
      fa = fb;
   }

   return out * xstep * 0.5;
}
vector ainteg(vector v[]; float xstep) {
   int    len = len(v), i;
   vector fa  = v[0];
   if(len<2) return fa;
   vector out = 0, fb;

   for(i=1; i<len; i++) {
      fb = v[i];
      out += fa + fb;
      fa = fb;
   }
   return out * xstep * 0.5;
}


//-------------------------------------------------------------------------------
// Spline Interpolation
//-------------------------------------------------------------------------------

// The bases for cubic splines, including the monotone-cubic and
// catmull-rom variants. They all assume row vectors (which is the order
// used in VEX).
//
// The hermite matrix can be used when working directly with tangents -- e.g:
// when the 4 "points" (paramters) are {p[i],p[i+1],m[i],m[i+1]}, with m[i] 
// and m[i+1] being the tangents associated with p[i] and p[i+1] respectively. 
// In that sense, it is the most general form as the rest of the cubic-spline 
// types really only represent different ways to auto-generate those tangents.
//
// | 1  t  t^2 t^3 |  |  1  0  0  0 | | p[i]   |
//                    |  0  0  1  0 | | p[i+1] |
//                    | -3  3 -2 -1 | | m[i]   |
//                    |  2 -2  1  1 | | m[i+1] |
//

#define SBASIS_HERMITE {  1,  0,  0,  0, \
                          0,  0,  1,  0, \
                         -3,  3, -2, -1, \
                          2, -2,  1,  1  }

// The catmull-rom spline is parameterized with {p[i-1],p[i],p[i+1],p[i+2]}
// and is identical to the Hermite above, except it chooses the tangents at
// the control points p[i] and p[i+1] to be:
//
//        p[i+1] - p[i-1]                  p[i+2] - p[i]
// m[i] = ---------------   and  m[i+1] =  -------------
//               2                               2
//
// which, when substituted into the Hermite form above, becomes:
//
// { p[i], p[i+1], (p[i+1]-p[i-1])/2, (p[i+2]-p[i])/2 } * SBASIS_HERMITE
//
// which can be written in matrix form as:
//
// | 1  t  t^2 t^3 |  |  1  0  0  0 | |  0    1   0   0   | | p[i-1] |
//                    |  0  0  1  0 | |  0    0   1   0   | | p[i]   |
//                    | -3  3 -2 -1 | | -1/2  0   1/2 0   | | p[i+1] |
//                    |  2 -2  1  1 | |  0   -1/2 0   1/2 | | p[i+2] |
//
// and pre-multiplying those two 4x4 matrices gives us the catrom basis:
//       |  0  2  0  0 |
//  1/2  | -1  0  1  0 |
//       |  2 -5  4 -1 |
//       | -1  3 -3  1 |
//

#define SBASIS_CATROM  {  0.0,  1.0,  0.0,  0.0, \
                         -0.5,  0.0,  0.5,  0.0, \
                          1.0, -2.5,  2.0, -0.5, \
                         -0.5,  1.5, -1.5,  0.5  }



// Interpolating functions supporting the set of spline types available in 
// standard ramps (including const and linear for completeness). 
// Naming convention: splerp_<segment_type>, to make it VOP friendly, where 
// the postfix (segment_type) can be resolved at preprocessing time. 
// The single API (branching on type) also included for calls from vex.

// constant
float  splerp_constant(float v0,v1, k) {
   return v0;
}
vector splerp_constantv(vector v0,v1; float k) {
   return v0;
}

// linear
float  splerp_linear(float v0,v1, k) {
   return lerp(v0,v1,k);
}
vector splerp_linearv(vector v0,v1; float k) {
   return lerp(v0,v1,k);
}

// Hermite
// These are here to support alternative tangent methods. E.g: to explore
// methods for achieving C1-continuity at non-uniform key intervals, which
// the methods below don't support (because they emulate VEX's spline()
// functions, which assume uniformly-spaced control points when solving
// cubic splines).
float  splerp_hermite(float p0,p1,m0,m1, k) {
   float k2 = k*k;
   vector4 w = set(1, k, k2, k2*k)*SBASIS_HERMITE;
   return w[0]*p0 + w[1]*p1 + w[2]*m0 + w[3]*m1;
}
vector splerp_hermitev(vector p0,p1,m0,m1; float k) {
   float k2 = k*k;
   vector4 w = set(1, k, k2, k2*k)*SBASIS_HERMITE;
   return w[0]*p0 + w[1]*p1 + w[2]*m0 + w[3]*m1;
}

// Catmull-Rom 
// These duplicate VEX's spline() behaviour (for catmull-rom segments), 
// which is *not* C1-continuous at the control points -- i.e: they don't 
// take into account the non-uniform spacing between control points.
float  splerp_catmullrom(float Pim1,Pi,Pip1,Pip2, k) {
   float k2 = k*k;
   vector4 w = set(1,k,k2,k2*k)*SBASIS_CATROM;
   return w[0]*Pim1 + w[1]*Pi + w[2]*Pip1 + w[3]*Pip2;
}
vector splerp_catmullromv(vector Pim1,Pi,Pip1,Pip2; float k) {
   float k2 = k*k;
   vector4 w = set(1,k,k2,k2*k)*SBASIS_CATROM;
   return w[0]*Pim1 + w[1]*Pi + w[2]*Pip1 + w[3]*Pip2;
}

// Monotone Cubic
// VEX's monotone-cubic interpolation is identical to splerp_hermite().
// The only difference is in the way the tangents are computed, but 
// the processing is the same (much like catmull-rom above is also hermite).
#define splerp_monotonecubic  splerp_hermite
#define splerp_monotonecubicv splerp_hermitev



//-------------------------------------------------------------------------------
// Slope calculation for monotone-cubic splines.
// Duplicates the monotone-cubic implementation in Houdini's ramps, which
// assumes evenly-spaced knots (!!!) -- i.e: these duplicate that limitation
// as well.
// They assume that the lookup position 'k' ('t' is reserved) is in [0,1]
//-------------------------------------------------------------------------------

// Monotone-Cubic: Calc the slope at the current knot 'p1', given the two 
// neighbouring knots 'p0' (previous knot) and 'p2' (next knot).
float ramp_slope_mc(float p0,p1,p2) 
{
   float m0  = p1-p0;
   float m1  = p2-p1;

   if(m0*m1<0 || equalzero(m0) || equalzero(m1))
      return 0;

   float am0 = abs(m0);
   float am1 = abs(m1);
   float l0  = am0 + 1;
   float l1  = am1 + 1;

   float w,slope;
   if(am1>=am0) {
      w     = (1.0 - (am0/am1)) / (1.0 + (l0/l1));
      slope = (1.0 + (2.0*w)) * m0;
   } else {
      w     = (1.0 - (am1/am0)) / (1.0 + (l1/l0));
      slope = (1.0 + (2.0*w)) * m1;
   }

   return slope;
}

// the vector version, slightly reworked -- should be a little faster than 
// calling the float verion three times...
vector ramp_slope_mc(vector p0,p1,p2) 
{
   vector m0  = p1-p0;
   vector m1  = p2-p1;

   vector am0 = abs(m0);
   vector am1 = abs(m1);
   vector l0  = am0 + 1;
   vector l1  = am1 + 1;

   float  w, m0i,m1i, a,b,c,d,e;
   vector slope;
   int    i;
   for(i=0;i<3;i++) 
   {
      m0i = m0[i];
      m1i = m1[i];
      if(m0i*m1i<0 || equalzero(m0i) || equalzero(m1i)) {
         slope[i] = 0;
      } else {
         if(am1[i]>=am0[i]) {
            a = am0[i]; b = am1[i]; c = l0[i]; d = l1[i]; e = m0i;
         } else {
            a = am1[i]; b = am0[i]; c = l1[i]; d = l0[i]; e = m1i;
         }
         w = (1.0 - (a/b)) / (1.0 + (c/d));
         slope[i] = (1.0 + (2.0*w)) * e;
      }
   }

   return slope;
}



//------------------------------------------------------------------------------
// Filter Width Calculation (per context)
//------------------------------------------------------------------------------
// For SOPs, POPs, and CHOPs
float  pyro_vopfw_SOP (float p)    { return 0; }
float  pyro_vopfw_SOP (vector2 p)   { return 0; }
float  pyro_vopfw_SOP (vector p)   { return 0; }
float  pyro_vopfw_SOP (vector4 p)  { return 0; }

float  pyro_vopfw_POP (float p)    { return 0; }
float  pyro_vopfw_POP (vector2 p)    { return 0; }
float  pyro_vopfw_POP (vector p)   { return 0; }
float  pyro_vopfw_POP (vector4 p)  { return 0; }

float  pyro_vopfw_CHOP (float p)   { return 0; }
float  pyro_vopfw_CHOP (vector2 p)  { return 0; }
float  pyro_vopfw_CHOP (vector p)  { return 0; }
float  pyro_vopfw_CHOP (vector4 p) { return 0; }

float  pyro_vopfw_VOP_CTXT (float p)   { return 0; }
float  pyro_vopfw_VOP_CTXT (vector2 p)  { return 0; }
float  pyro_vopfw_VOP_CTXT (vector p)  { return 0; }
float  pyro_vopfw_VOP_CTXT (vector4 p) { return 0; }

// For the shading contexts and COPs
float pyro_vopfw_SHADING(float p) {
   return max(max(abs(Du(p)),abs(Dv(p))),1e-6);
   //return max(max(abs(Du(p)),max(abs(Dv(p)),abs(Dw(p)))),1e-6);
   //return abs(Du(p))+abs(Dv(p))+abs(Dw(p));
   //return max(length(set(abs(Du(p)),abs(Dv(p)),abs(Dw(p)))),1e-6);
}
#define pyro_vopfw_COP  pyro_vopfw_SHADING
#define pyro_vopfw_COP2 pyro_vopfw_SHADING

float pyro_vopfw_SHADING(vector2 p) {
   //return max(max(length(Du(p)),length(Dv(p))),1e-6);
   return max(max(abs(Du(vector(p)))+abs(Dv((vector)p))),1e-6);
}
#define pyro_vopfwv_COP  pyro_vopfw_SHADING
#define pyro_vopfwv_COP2 pyro_vopfw_SHADING

float pyro_vopfw_SHADING(vector p) {
   //return max(max(length(Du(p)),length(Dv(p))),1e-6);
   return max(max(abs(Du(p))+abs(Dv(p))),1e-6);
}
#define pyro_vopfwv_COP  pyro_vopfw_SHADING
#define pyro_vopfwv_COP2 pyro_vopfw_SHADING

float pyro_vopfw_SHADING(vector4 p) {
   return max(max(length(Du(p)),length(Dv(p))),1e-6);
}
#define pyro_vopfw_COP  pyro_vopfw_SHADING
#define pyro_vopfw_COP2 pyro_vopfw_SHADING

// The front-end API to all the above vop filter functions
#define VOPFW(p)  CAT3(pyro_vopfw,_,VOP_CTXT)(p)

//------------------------------------------------------------------------------
// Contour and Softclip
//------------------------------------------------------------------------------

// Soft Clip - reduces a field's high-intensity range by
// smoothly compressing it into a log scale (beyond a given threshold).
// A smoother version of a "ceiling clamp".
float softclip (
      float f;       // field value
      float start;   // start val for compression
      float c;       // compression
   )
{
   float out = f;
   if(f>start && c>0) {
      float ki = 1.0 / c;
      float w  = 1.0 / (c*log(10.0));
      float k  = log10(pow(w,ki));
      out = log10(pow((f-start)+w,ki)) - k + start;
   }
   return out;
}

vector softclip(vector x; float start, c) {
   return set(softclip(x[0],start,c),
              softclip(x[1],start,c),
              softclip(x[2],start,c));
}
vector4 softclip(vector4 x; float start, c) {
   return set(softclip(x[0],start,c),
              softclip(x[1],start,c),
              softclip(x[2],start,c),
              softclip(x[3],start,c));
}


// Contour - adds contrast to a field's low-intensity range by
// smoothly attenuating values below a given threshold.
// A smoother version of a "floor clamp".
// Implemented as 3x^(2*p) - 2x^(3*p), where 'p' is the 'rolloff'.
// NOTE: This is not quite the same as VEX's smooth() functions where 'p'
// (the 'rolloff') is implemented as (3x^2 - 2x^3)^p, which gives
// a slightly different rate of "rolling" when p!=1, but this version
// is otherwise identical, and more importantly, much easier to 
// integrate!. They are indistiguishable from the point of view of
// usability, but *not* identical for p!=1.


// contour, unfiltered
float contour (
      float f;       // field value
      float e;       // field value interpreted as the "contour" edge
      float sharp;   // sharpening amount [0,+inf] (practical range [0,100])
   )
{
   if(f<0) return 0;
   if(f>e || sharp<=0) return f;
   float ee = max(1e-3,abs(e));
   float x = f/ee;
   if(sharp==1) return f * (3*x*x - 2*x*x*x);
   return f * (3*pow(x,2*sharp) - 2*pow(x,3*sharp));
}
vector contour(vector f; float e, sharp) {
   return set(contour(f.x,e,sharp),
              contour(f.y,e,sharp),
              contour(f.z,e,sharp) );
}
vector4 contour(vector4 f; float e, sharp) {
   return set(contour(f.x,e,sharp),
              contour(f.y,e,sharp),
              contour(f.z,e,sharp),
              contour(f.w,e,sharp) );
}


// Integrates the contour() function [= 3x^(2*p) - 2x^(3*p)] in the
// interval [a,b] where a and b are in [0,1] and b>=a.
// Parts of the interval [a,b] that fall outside [0,1] are trivially 0,
// and b-1 respectively, and they're dealt with by the caller.
float contourI(float a,b,p) {
   if(p==1) return a*a*a*(a*0.5-1) + b*b*b*(1-b*0.5);
   float m = 1 + 2*p;
   float n = 1 + 3*p;
   float am = pow(a,m), an = pow(a,n);
   float bm = pow(b,m), bn = pow(b,n);
   float out = -3*am + 2*an + 3*bm - 2*bn - 9*am*p + 4*an*p +9*bm*p - 4*bn*p;
   out /= m*n;
   return max(0,out);
}

// filtered contour
float contour(float f, e, sharp, fw) {
   if(e<=0) return 0;
   if(fw<C_FLOATEPSILON) return contour(f,e,sharp);

   float ee = max(1e-3,abs(e));
   float w  = 0;
   float x0 = f-(fw*.5);
   float x1 = x0+fw;
   if(x0>=e) {
      w = 1;
   } else if(x1<0) {
      w = 0;
   } else {
      float ss = 1/ee;
      x0 *= ss;
      x1 *= ss;
      if(x0<1 && x1>0) w += contourI(max(0,x0),min(1,x1),sharp);
      if(x1>1) w += x1-1; 
      w /= ss*fw;
   }
   return f * w;
}
vector contour(vector f; float e, sharp, fw) {
   return set(contour(f.x,e,sharp,fw),
              contour(f.y,e,sharp,fw),
              contour(f.z,e,sharp,fw) );
}
vector4 contour(vector4 f; float e, sharp, fw) {
   return set(contour(f.x,e,sharp,fw),
              contour(f.y,e,sharp,fw),
              contour(f.z,e,sharp,fw),
              contour(f.w,e,sharp,fw) );
}

// Compute the mask for a field.  For micropolygon renders, this will
// enlarge by the filter width.
float pyro_fieldmask(float val; float edge)
{
    float fw = edge;
    // TODO: This introduces derivatives of density into default pyro
    // renders, even if it is ray traced, leading to substantial
    // performance degredation due to the increased number of field
    // evaluations.
#if 0
    if (!israytracing())
	fw = VOPFW(val);
#endif
    return filterstep(edge, max(val-fw, 0), val,
	    "filter", "gauss", "width", 2);
}


//-------------------------------------------------------------------------------
// Field Functions
//-------------------------------------------------------------------------------
vector fieldgradient(vector dfield; vector dspace[]) {
   return dfield[0]*dspace[0] + 
          dfield[1]*dspace[1] + 
          dfield[2]*dspace[2];
}

vector fieldgradientv(vector dfield[]; vector dspace[]) {
   return dot(dfield[0],dspace[0])*dspace[0] + 
          dot(dfield[1],dspace[1])*dspace[1] + 
          dot(dfield[2],dspace[2])*dspace[2];
}




//-------------------------------------------------------------------------------
// Assorted VOP Utilities
//-------------------------------------------------------------------------------

// Re-Interpolation of rest ratios for the dualrest vop
#define RRINTERP_LINEAR    0
#define RRINTERP_CUBIC     1
#define RRINTERP_COSINE    2
float rrinterp(float lin; int newinterp) {
   if(newinterp==RRINTERP_CUBIC)  return smooth(0,1,lin);
   if(newinterp==RRINTERP_COSINE) return 1 - (cos(lin*C_PI)*0.5+0.5);
   return lin;
}

#endif // End pyro_utils_h
