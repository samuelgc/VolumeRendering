/****************************************************************************
   Description:  Support for filtered ramps.
                 Also extends standard ramps to support various boundary
                 conditions. All combinations split into different functions
                 for VOP calls (no branching).

*****************************************************************************/

#ifndef pyro_aaramp_h__GUARD
#define pyro_aaramp_h__GUARD


#include <pyro_math.h>
#include <pyro_color.h>
#include <pyro_utils.h>



// Ramp Condition
// Ensures the given ramp has keys at both limits {0,1}. Returns the 
// arrays of the conditioned ramp.
// This helps reduce some of the branching during analysis -> faster compiles.
//-------------------------------------------------------------------------------
void ramp_condition(
      string inbases[];        // input ramp: basis strings
      float  inkeys[];         // input ramp: key positions
      float  invals[];         // input ramp: key values
      export string bases[];   // output ramp: basis strings
      export float  keys[];    // output ramp: key positions
      export float  vals[];    // output ramp: key values
   )
{
   int ink      = len(inkeys);
   int inl      = max(0,ink-1);
   int addfirst = (inkeys[0]>0);
   int addlast  = (inkeys[inl]<1);

   if(!(addfirst|addlast)) {
      bases = inbases;
      keys  = inkeys;
      vals  = invals;
      return;
   }

   int onk = ink+addfirst+addlast;
   resize(bases,onk);
   resize(keys,onk);
   resize(vals,onk);

   int k = 0;
   if(addfirst) {
      keys[0]  = 0;
      bases[0] = "constant";
      vals[0]  = invals[0];
      ++k;
   }

   int i; 
   for(i=0;i<ink;i++) {
      keys[k]  = inkeys[i];
      bases[k] = inbases[i];
      vals[k]  = invals[i];
      ++k;
   }

   if(addlast) {
      bases[k-1] = "constant";
      keys[k]    = 1;
      bases[k]   = "constant";
      vals[k]    = invals[inl];
   }

}
void ramp_condition(
      string inbases[];        // input ramp: basis strings
      float  inkeys[];         // input ramp: key positions
      vector invals[];         // input ramp: key values
      export string bases[];   // output ramp: basis strings
      export float  keys[];    // output ramp: key positions
      export vector vals[];    // output ramp: key values
   )
{
   int ink      = len(inkeys);
   int inl      = max(0,ink-1);
   int addfirst = (inkeys[0]>0);
   int addlast  = (inkeys[inl]<1);

   if(!(addfirst|addlast)) {
      bases = inbases;
      keys  = inkeys;
      vals  = invals;
      return;
   }

   int onk = ink+addfirst+addlast;
   resize(bases,onk);
   resize(keys,onk);
   resize(vals,onk);

   int k = 0;
   if(addfirst) {
      keys[0]  = 0;
      bases[0] = "constant";
      vals[0]  = invals[0];
      ++k;
   }

   int i; 
   for(i=0;i<ink;i++) {
      keys[k]  = inkeys[i];
      bases[k] = inbases[i];
      vals[k]  = invals[i];
      ++k;
   }

   if(addlast) {
      bases[k-1] = "constant";
      keys[k]    = 1;
      bases[k]   = "constant";
      vals[k]    = invals[inl];
   }

}

// Transform to Linear-RGB space
// For color ramps, both the color model (XYZ, RGB, HSL, etc) and the
// gamma are inspected. For float ramps, only the gamma setting is applied.
// In both cases, the transform is applide in-place.
//-------------------------------------------------------------------------------
vector ramp_torgb(string space; vector col) {
   if(space!="rgb") {
      if(space=="XYZ") return xyztorgbgamut(col);
      return ctransform("cspace:"+space,"cspace:rgb",col);
   }
   return col;
}
float ramp_tolinear(string space; float val; float customgamma) {
   if(space=="linear")  return val;
   if(space=="srgb")    return fromSRGB(val);
   float sg = sign(customgamma);
   return pow(val,max(1e-6,customgamma*sg)*sg);
}
vector ramp_tolinear(string space; vector val; float customgamma) {
   if(space=="linear")  return val;
   if(space=="srgb")    return fromSRGB(val);
   float sg = sign(customgamma);
   return pow(val,max(1e-6,customgamma*sg)*sg);
}
vector4 ramp_tolinear(string space; vector4 val; float customgamma) {
   if(space=="linear")  return val;
   if(space=="srgb")    return fromSRGB(val);
   float sg = sign(customgamma);
   return pow(val,max(1e-6,customgamma*sg)*sg);
}

vector ramp_rgbtorgb(vector c) { return c; }    // safety

vector ramp_hsvtorgb(vector c) {
   return ctransform("cspace:hsv","cspace:rgb",c);
}

vector ramp_labtorgb(vector c) { 
   return ctransform("cspace:Lab","cspace:rgb",c);
}

vector ramp_hsltorgb(vector c) { 
   return ctransform("cspace:hsl","cspace:rgb",c);
}

vector ramp_xyztorgb(vector c) {
   //return ctransform("cspace:XYZ","cspace:rgb",c);
   return xyztorgbgamut(c);
}

vector ramp_fromlinear(vector c; float g) { return c; }  // safety
float  ramp_fromlinear(float  c; float g) { return c; }  // safety

vector ramp_fromsrgb(vector c; float g) {
   return fromSRGB(c);
}
float ramp_fromsrgb(float c; float g) {
   return fromSRGB(c);
}
vector4 ramp_fromsrgb(vector4 c; float g) {
   return fromSRGB(c);
}

vector ramp_fromcustom(vector c; float g) {
   return pow(c,g);
}
float ramp_fromcustom(float c; float g) {
   return pow(c,g);
}
vector4 ramp_fromcustom(vector4 c; float g) {
   return pow(c,g);
}


//-------------------------------------------------------------------------------
// Ramp Lookup
// Returns a ramp's value at lookup position 'k' in [-inf,+inf];
// Implements boundary conditions for left/right boundaries independently.
// Boundary types are strings and can be one of:
// "hold"   -  clamp boundary value
// "cycle"  -  repeat ramp
// "accum"  -  repeat and accumulate (first val shifted to start at last val)
// "mirror" -  cycle while alternating directions
// "slope"  -  continue along function slope at boundary
//-------------------------------------------------------------------------------

// delta step for calculating slope
#define  MDELTA (1e-3)
#define IMDELTA (1e+3)

// These custom versions are here for debugging purposes only.
// Their output should match those of the native VEX versions (below), as
// these custom algorithms are used when computing a ramp's integral for AA.
float ramp_lookup_debug(float k; string bases[]; float keys[]; float vals[])
{
   int n  = len(keys);
   if(n<2) return vals[0];
   int nl = n-1;

   float f  = spline("solvelinear",k,keys) * nl;
   int   n1 = floor(f);
   f -= n1;

   if(n1==nl) return vals[nl];
   if(bases[n1]=="constant") 
      return vals[n1];

   int   n2 = min(nl,n1+1);
   if(bases[n1]=="linear") 
      return lerp(vals[n1],vals[n2],f);

   int   n0 = max(0,n1-1);
   int   n3 = min(nl,n2+1);
   float p0 = vals[n0], p1 = vals[n1], p2 = vals[n2], p3 = vals[n3];

   if(bases[n1]=="catmull-rom")
      return splerp_catmullrom(p0,p1,p2,p3,f);

   float m1 = ramp_slope_mc(p0,p1,p2);
   float m2 = ramp_slope_mc(p1,p2,p3);
   return splerp_hermite(p1,p2,m1,m2,f);
}
vector ramp_lookup_debug(float k; string bases[]; float keys[]; vector vals[])
{
   int n  = len(keys);
   if(n<2) return vals[0];
   int nl = n-1;

   float  f  = spline("solvelinear",k,keys) * nl;
   int    n1 = floor(f);
   f -= n1;

   if(n1==nl) return vals[nl];
   if(bases[n1]=="constant") 
      return vals[n1];

   int    n2 = min(nl,n1+1);
   if(bases[n1]=="linear") 
      return lerp(vals[n1],vals[n2],f);

   int    n0 = max(0,n1-1);
   int    n3 = min(nl,n2+1);
   vector p0 = vals[n0], p1 = vals[n1], p2 = vals[n2], p3 = vals[n3];

   if(bases[n1]=="catmull-rom")
      return splerp_catmullromv(p0,p1,p2,p3,f);

   float  spacing = 1.0 / (float)nl;
   vector m1 = ramp_slope_mc(p0,p1,p2);
   vector m2 = ramp_slope_mc(p1,p2,p3);
   return splerp_hermitev(p1,p2,m1,m2,f);
}

// These are the native versions which are actually called to perform a
// point sample (single value lookup). They assume k is in [0,1]
float ramp_lookup(float k; string bases[]; float keys[]; float vals[])
{
   return spline(bases,spline("solvelinear",k,keys),vals);
}
vector ramp_lookup(float k; string bases[]; float keys[]; vector vals[])
{
   return spline(bases,spline("solvelinear",k,keys),vals);
}



//-------------------------------------------------------------------------------
// Determine the slope at either end of a ramp
//-------------------------------------------------------------------------------

// differential step (dx) for calculating slope
#define  MDELTA (1e-3)
#define IMDELTA (1e+3)


void ramp_slopes (
      int    nl;
      string bases[]; 
      float  keys[]; 
      float  vals[];
      export float lslope;
      export float rslope;
   )
{
   lslope = (ramp_lookup(MDELTA,bases,keys,vals) - vals[0]) * IMDELTA;
   rslope = (vals[nl] - ramp_lookup(1-MDELTA,bases,keys,vals)) * IMDELTA;
}
void ramp_slopes (
      int    nl;
      string bases[]; 
      float  keys[]; 
      vector vals[];
      export vector lslope;
      export vector rslope;
   )
{
   lslope = (ramp_lookup(MDELTA,bases,keys,vals) - vals[0]) * IMDELTA;
   rslope = (vals[nl] - ramp_lookup(1-MDELTA,bases,keys,vals)) * IMDELTA;
}



//-------------------------------------------------------------------------------
// Ramp Filters
//-------------------------------------------------------------------------------

#define RAMP_FILTERED_PARMS(T) \
   int nl; float x,a,b; T lv,rv,lm,rm,dv; \
   string bases[]; float keys[]; T vals[];

#define RAMP_FILTERED_ARGS \
   nl, x,a,b, lv,rv,lm,rm,dv, \
   bases, keys, vals

float mirrorindex(float x) {
   return 1.0 - 2.0*abs(0.5-modulo(x*.5,1));
}
vector mirrorlimits(int seg; float a, b) {
   return (seg&1) ? set(1-b,1-a,0) : set(a,b,0);
}

// FILTER: Point
//--------------------------------------
// Float
float ramp_hold_hold_point( RAMP_FILTERED_PARMS(float) ) {
   return ramp_lookup(x,bases,keys,vals);
}
float ramp_hold_slope_point( RAMP_FILTERED_PARMS(float) ) {
   if(x<=1) return ramp_lookup(x,bases,keys,vals);
   return rv + (x-1)*rm;
}
float ramp_hold_cycle_point( RAMP_FILTERED_PARMS(float) ) {
   if(x<=1) return ramp_lookup(x,bases,keys,vals);
   return ramp_lookup(modulo(x,1),bases,keys,vals);
}
float ramp_hold_mirror_point( RAMP_FILTERED_PARMS(float) ) {
   if(x<=1) return ramp_lookup(x,bases,keys,vals);
   return ramp_lookup(mirrorindex(x),bases,keys,vals);
}
float ramp_hold_accum_point( RAMP_FILTERED_PARMS(float) ) {
   if(x<=1) return ramp_lookup(x,bases,keys,vals);
   return ramp_lookup(modulo(x,1),bases,keys,vals) + dv*floor(x);
}

float ramp_slope_hold_point( RAMP_FILTERED_PARMS(float) ) {
   if(x<0) return lv + x*lm;
   return ramp_lookup(x,bases,keys,vals);
}
float ramp_slope_slope_point( RAMP_FILTERED_PARMS(float) ) {
   if(x<0) return lv + x*lm;
   if(x>1) return rv + (x-1)*rm;
   return ramp_lookup(x,bases,keys,vals);
}
float ramp_slope_cycle_point( RAMP_FILTERED_PARMS(float) ) {
   if(x<0) return lv + x*lm;
   return ramp_lookup(modulo(x,1),bases,keys,vals);
}
float ramp_slope_mirror_point( RAMP_FILTERED_PARMS(float) ) {
   if(x<0) return lv + x*lm;
   return ramp_lookup(mirrorindex(x),bases,keys,vals);
}
float ramp_slope_accum_point( RAMP_FILTERED_PARMS(float) ) {
   if(x<0) return lv + x*lm;
   return ramp_lookup(modulo(x,1),bases,keys,vals) + dv*floor(x);
}

float ramp_cycle_hold_point( RAMP_FILTERED_PARMS(float) ) {
   if(x<=1) return ramp_lookup(modulo(x,1),bases,keys,vals);
   return ramp_lookup(x,bases,keys,vals);
}
float ramp_cycle_slope_point( RAMP_FILTERED_PARMS(float) ) {
   if(x<=1) return ramp_lookup(modulo(x,1),bases,keys,vals);
   return rv + (x-1)*rm;
}
float ramp_cycle_cycle_point( RAMP_FILTERED_PARMS(float) ) {
   return ramp_lookup(modulo(x,1),bases,keys,vals);
}
float ramp_cycle_mirror_point( RAMP_FILTERED_PARMS(float) ) {
   if(x<=1) return ramp_lookup(modulo(x,1),bases,keys,vals);
   return ramp_lookup(mirrorindex(x),bases,keys,vals);
}
float ramp_cycle_accum_point( RAMP_FILTERED_PARMS(float) ) {
   if(x<=1) return ramp_lookup(modulo(x,1),bases,keys,vals);
   return ramp_lookup(modulo(x,1),bases,keys,vals) + dv*floor(x);
}

float ramp_mirror_hold_point( RAMP_FILTERED_PARMS(float) ) {
   if(x<=1) return ramp_lookup(mirrorindex(x),bases,keys,vals);
   return rv;
}
float ramp_mirror_slope_point( RAMP_FILTERED_PARMS(float) ) {
   if(x<=1) return ramp_lookup(mirrorindex(x),bases,keys,vals);
   return rv + (x-1)*rm;
}
float ramp_mirror_cycle_point( RAMP_FILTERED_PARMS(float) ) {
   if(x<0) return ramp_lookup(mirrorindex(x),bases,keys,vals);
   return ramp_lookup(modulo(x,1),bases,keys,vals);
}
float ramp_mirror_mirror_point( RAMP_FILTERED_PARMS(float) ) {
   return ramp_lookup(mirrorindex(x),bases,keys,vals);
}
float ramp_mirror_accum_point( RAMP_FILTERED_PARMS(float) ) {
   if(x<0) return ramp_lookup(mirrorindex(x),bases,keys,vals);
   return ramp_lookup(modulo(x,1),bases,keys,vals) + dv*floor(x);
}

float ramp_accum_hold_point( RAMP_FILTERED_PARMS(float) ) {
   if(x<0) return ramp_lookup(modulo(x,1),bases,keys,vals) + dv*floor(x);
   return ramp_lookup(x,bases,keys,vals);
}
float ramp_accum_slope_point( RAMP_FILTERED_PARMS(float) ) {
   if(x<=1) return ramp_lookup(modulo(x,1),bases,keys,vals) + dv*floor(x);
   return rv + (x-1)*rm;
}
float ramp_accum_cycle_point( RAMP_FILTERED_PARMS(float) ) {
   if(x<=1) return ramp_lookup(modulo(x,1),bases,keys,vals) + dv*floor(x);
   return ramp_lookup(modulo(x,1),bases,keys,vals);
}
float ramp_accum_mirror_point( RAMP_FILTERED_PARMS(float) ) {
   if(x<=1) return ramp_lookup(modulo(x,1),bases,keys,vals) + dv*floor(x);
   return ramp_lookup(mirrorindex(x),bases,keys,vals);
}
float ramp_accum_accum_point( RAMP_FILTERED_PARMS(float) ) {
   return ramp_lookup(modulo(x,1),bases,keys,vals) + dv*floor(x);
}


// Vector
vector ramp_hold_hold_point( RAMP_FILTERED_PARMS(vector) ) {
   return ramp_lookup(x,bases,keys,vals);
}
vector ramp_hold_slope_point( RAMP_FILTERED_PARMS(vector) ) {
   if(x<=1) return ramp_lookup(x,bases,keys,vals);
   return rv + (x-1)*rm;
}
vector ramp_hold_cycle_point( RAMP_FILTERED_PARMS(vector) ) {
   if(x<=1) return ramp_lookup(x,bases,keys,vals);
   return ramp_lookup(modulo(x,1),bases,keys,vals);
}
vector ramp_hold_mirror_point( RAMP_FILTERED_PARMS(vector) ) {
   if(x<=1) return ramp_lookup(x,bases,keys,vals);
   return ramp_lookup(mirrorindex(x),bases,keys,vals);
}
vector ramp_hold_accum_point( RAMP_FILTERED_PARMS(vector) ) {
   if(x<=1) return ramp_lookup(x,bases,keys,vals);
   return ramp_lookup(modulo(x,1),bases,keys,vals) + dv*floor(x);
}

vector ramp_slope_hold_point( RAMP_FILTERED_PARMS(vector) ) {
   if(x<0) return lv + x*lm;
   return ramp_lookup(x,bases,keys,vals);
}
vector ramp_slope_slope_point( RAMP_FILTERED_PARMS(vector) ) {
   if(x<0) return lv + x*lm;
   if(x>1) return rv + (x-1)*rm;
   return ramp_lookup(x,bases,keys,vals);
}
vector ramp_slope_cycle_point( RAMP_FILTERED_PARMS(vector) ) {
   if(x<0) return lv + x*lm;
   return ramp_lookup(modulo(x,1),bases,keys,vals);
}
vector ramp_slope_mirror_point( RAMP_FILTERED_PARMS(vector) ) {
   if(x<0) return lv + x*lm;
   return ramp_lookup(mirrorindex(x),bases,keys,vals);
}
vector ramp_slope_accum_point( RAMP_FILTERED_PARMS(vector) ) {
   if(x<0) return lv + x*lm;
   return ramp_lookup(modulo(x,1),bases,keys,vals) + dv*floor(x);
}

vector ramp_cycle_hold_point( RAMP_FILTERED_PARMS(vector) ) {
   if(x<=1) return ramp_lookup(modulo(x,1),bases,keys,vals);
   return ramp_lookup(x,bases,keys,vals);
}
vector ramp_cycle_slope_point( RAMP_FILTERED_PARMS(vector) ) {
   if(x<=1) return ramp_lookup(modulo(x,1),bases,keys,vals);
   return rv + (x-1)*rm;
}
vector ramp_cycle_cycle_point( RAMP_FILTERED_PARMS(vector) ) {
   return ramp_lookup(modulo(x,1),bases,keys,vals);
}
vector ramp_cycle_mirror_point( RAMP_FILTERED_PARMS(vector) ) {
   if(x<=1) return ramp_lookup(modulo(x,1),bases,keys,vals);
   return ramp_lookup(mirrorindex(x),bases,keys,vals);
}
vector ramp_cycle_accum_point( RAMP_FILTERED_PARMS(vector) ) {
   if(x<=1) return ramp_lookup(modulo(x,1),bases,keys,vals);
   return ramp_lookup(modulo(x,1),bases,keys,vals) + dv*floor(x);
}

vector ramp_mirror_hold_point( RAMP_FILTERED_PARMS(vector) ) {
   if(x<=1) return ramp_lookup(mirrorindex(x),bases,keys,vals);
   return rv;
}
vector ramp_mirror_slope_point( RAMP_FILTERED_PARMS(vector) ) {
   if(x<=1) return ramp_lookup(mirrorindex(x),bases,keys,vals);
   return rv + (x-1)*rm;
}
vector ramp_mirror_cycle_point( RAMP_FILTERED_PARMS(vector) ) {
   if(x<0) return ramp_lookup(mirrorindex(x),bases,keys,vals);
   return ramp_lookup(modulo(x,1),bases,keys,vals);
}
vector ramp_mirror_mirror_point( RAMP_FILTERED_PARMS(vector) ) {
   return ramp_lookup(mirrorindex(x),bases,keys,vals);
}
vector ramp_mirror_accum_point( RAMP_FILTERED_PARMS(vector) ) {
   if(x<0) return ramp_lookup(mirrorindex(x),bases,keys,vals);
   return ramp_lookup(modulo(x,1),bases,keys,vals) + dv*floor(x);
}

vector ramp_accum_hold_point( RAMP_FILTERED_PARMS(vector) ) {
   if(x<0) return ramp_lookup(modulo(x,1),bases,keys,vals) + dv*floor(x);
   return ramp_lookup(x,bases,keys,vals);
}
vector ramp_accum_slope_point( RAMP_FILTERED_PARMS(vector) ) {
   if(x<=1) return ramp_lookup(modulo(x,1),bases,keys,vals) + dv*floor(x);
   return rv + (x-1)*rm;
}
vector ramp_accum_cycle_point( RAMP_FILTERED_PARMS(vector) ) {
   if(x<=1) return ramp_lookup(modulo(x,1),bases,keys,vals) + dv*floor(x);
   return ramp_lookup(modulo(x,1),bases,keys,vals);
}
vector ramp_accum_mirror_point( RAMP_FILTERED_PARMS(vector) ) {
   if(x<=1) return ramp_lookup(modulo(x,1),bases,keys,vals) + dv*floor(x);
   return ramp_lookup(mirrorindex(x),bases,keys,vals);
}
vector ramp_accum_accum_point( RAMP_FILTERED_PARMS(vector) ) {
   return ramp_lookup(modulo(x,1),bases,keys,vals) + dv*floor(x);
}








// FILTER: Box
//--------------------------------------

// native width of this filter
float fscale_box() { return 1; }

// convolve with a constant -- for "hold" boundaries
float ramp_sumconst_box(float val,size) {
   return val*size;
}
vector ramp_sumconst_box(vector val; float size) {
   return val*size;
}

// convolve with a line -- for "slope" boundaries
float ramp_sumslope_box(float val,slope,a,b) {
   return 0.5 * (b-a) * (2*val + slope*(a+b));
}
vector ramp_sumslope_box(vector val,slope; float a,b) {
   return 0.5 * (b-a) * (2*val + slope*(a+b));
}

// Convolve box filter with a portion {a,b} of a single ramp segment (at index n)
float ramp_sumseg_box(
      int    n;         // array index at segment start
      int    nl;        // last array index
      float  ci;        // const of integration
      float  a,b;       // limits of integration ({a,b} in [0,1], a<=b)
      string bases[];   // ramp basis strings
      float  keys[];    // ramp key positions
      float  vals[]     // ramp key values
   ) 
{
   int   n0  = max(0,n-1), n1 = n, n2 = min(nl,n+1), n3 = min(nl,n2+1);
   float len = keys[n2]-keys[n1];
   float v0  = vals[n0]+ci, v1 = vals[n1]+ci, 
         v2  = vals[n2]+ci, v3 = vals[n3]+ci;

   if(bases[n]=="constant")
      return len * (b - a) * v1;

   if(bases[n]=="linear")
      return len * 0.5*(a - b)*( (-2 + a + b)*v1 - (a + b)*v2 );

   // catmull-rom (in Horner form)
   if(bases[n]=="catmull-rom") {
      return (len/24) * ( 
         a*(-24*v1+a*(6*v0-6*v2+a*(-8*v0+20*v1-16*v2+a*(3*v0-9*v1+9*v2-3*v3)+4*v3)))
        +b*(24*v1+b*(-6*v0+6*v2+b*(8*v0-20*v1+16*v2-4*v3+b*(-3*v0+9*v1-9*v2+3*v3)))) );
   }

   // monotone cubic (in Horner form)
   float m1 = ramp_slope_mc(v0,v1,v2);
   float m2 = ramp_slope_mc(v1,v2,v3);
   return (len/12) * (
        a*(-12*v1+a*(-6*m1+a*(8*m1+4*m2+12*v1-12*v2+a*(-3*m1-3*m2-6*v1+6*v2))))
       +b*(12*v1+b*(6*m1+b*(-8*m1-4*m2-12*v1+b*(3*m1+3*m2+6*v1-6*v2)+12*v2))) );

}

vector ramp_sumseg_box(
      int    n;         // array index at segment start
      int    nl;        // last array index
      vector ci;        // const of integration
      float  a,b;       // limits of integration ({a,b} in [0,1], a<=b)
      string bases[];   // ramp basis strings
      float  keys[];    // ramp key positions
      vector vals[]     // ramp key values
   ) 
{
   int    n0 = max(0,n-1), n1 = n, n2 = min(nl,n+1), n3 = min(nl,n2+1);
   float  len = keys[n2]-keys[n1];
   vector v0 = vals[n0]+ci, v1 = vals[n1]+ci, 
          v2 = vals[n2]+ci, v3 = vals[n3]+ci;

   if(bases[n]=="constant") {
      return len * (b - a) * v1;
   }

   if(bases[n]=="linear")
      return len * 0.5*(a - b)*( (-2 + a + b)*v1 - (a + b)*v2 );

   // catmull-rom (in Horner form)
   if(bases[n]=="catmull-rom") {
      return (len/24) * ( 
         a*(-24*v1+a*(6*v0-6*v2+a*(-8*v0+20*v1-16*v2+a*(3*v0-9*v1+9*v2-3*v3)+4*v3)))
        +b*(24*v1+b*(-6*v0+6*v2+b*(8*v0-20*v1+16*v2-4*v3+b*(-3*v0+9*v1-9*v2+3*v3)))) );
   }

   // monotone cubic (in Horner form)
   vector m1 = ramp_slope_mc(v0,v1,v2);
   vector m2 = ramp_slope_mc(v1,v2,v3);
   return (len/12) * (
        a*(-12*v1+a*(-6*m1+a*(8*m1+4*m2+12*v1-12*v2+a*(-3*m1-3*m2-6*v1+6*v2))))
       +b*(12*v1+b*(6*m1+b*(-8*m1-4*m2-12*v1+b*(3*m1+3*m2+6*v1-6*v2)+12*v2))) );

}


// convolve box filter with a portion {a,b} of the entire [0,1] ramp
float ramp_sum_box (
      int    nl;        // last array index
      float  ci;        // const of integration
      float  a,b;       // limits of integration in [0,1]
      string bases[];   // ramp basis strings
      float  keys[];    // ramp key positions
      float  vals[];    // ramp key values
   ) 
{
   float sum = 0;

   float f0 = (spline("solvelinear",a,keys));
   float f1 = (spline("solvelinear",b,keys));
   int   n0 = floor(f0*nl);
   int   n1 = (int)ceil(f1*nl);
   int n; float k0,k1,il,sa,sb;
   for(n=n0;n<n1;n++) { 
      k0 = keys[n];
      k1 = keys[n+1];
      il = 1/(k1-k0);
      sa = (max(k0,a)-k0)*il;  // segment limits {sa,sb} in [0,1]
      sb = (min(k1,b)-k0)*il;
      sum += ramp_sumseg_box(n,nl, ci, sa,sb, bases,keys,vals);
   }

   return sum;
}

vector ramp_sum_box (
      int    nl;        // last array index
      vector ci;        // const of integration
      float  a,b;       // limits of integration in [0,1]
      string bases[];   // ramp basis strings
      float  keys[];    // ramp key positions
      vector vals[];    // ramp key values
   ) 
{
   vector sum = 0;

   float f0 = (spline("solvelinear",a,keys));
   float f1 = (spline("solvelinear",b,keys));
   int   n0 = floor(f0*nl);
   int   n1 = (int)ceil(f1*nl);
   int n; float k0,k1,il,sa,sb;
   for(n=n0;n<n1;n++) { 
      k0 = keys[n];
      k1 = keys[n+1];
      il = 1/(k1-k0);
      sa = (max(k0,a)-k0)*il;  // segment limits {sa,sb} in [0,1]
      sb = (min(k1,b)-k0)*il;
      sum += ramp_sumseg_box(n,nl, ci, sa,sb, bases,keys,vals);
   }

   return sum;
}

#define SUM_BOX_PARMS(T) \
   int nl; float a,b; T dv; string bases[]; float keys[]; T vals[];
#define SUM_BOX_ARGS(T,a,b) \
   nl, a,b, dv, bases, keys, vals

// Cycle
float  ramp_sumcycle_box( SUM_BOX_PARMS(float ) ) {
   float sum = 0;
   float fa = floor(a), ca = fa+1, fb = floor(b);
   // first possibly partial cell
   sum += ramp_sum_box(nl,0,a-fa,min(1,b-fa),bases,keys,vals);
   if(fb>fa) {
      // last possibly partial cell
      sum += ramp_sum_box(nl,0,0,b-fb,bases,keys,vals);
      // plus all the complete cells in between (there are 'nr' of them)
      int nr = (int)(fb-ca);
      if(nr>0) sum += nr * ramp_sum_box(nl,0,0,1,bases,keys,vals);
   }
   return sum;
}

vector ramp_sumcycle_box( SUM_BOX_PARMS(vector) ) {
   vector sum = 0;
   float fa = floor(a), ca = fa+1, fb = floor(b);
   sum += ramp_sum_box(nl,VZERO,a-fa,min(1,b-fa),bases,keys,vals);
   if(fb>fa) {
      sum += ramp_sum_box(nl,VZERO,0,b-fb,bases,keys,vals);
      int nr = (int)(fb-ca);
      if(nr>0) sum += nr * ramp_sum_box(nl,VZERO,0,1,bases,keys,vals);
   }
   return sum;
}

// Mirror
float  ramp_summirror_box( SUM_BOX_PARMS(float ) ) {
   float  sum = 0;
   // Same as cycle, except limits are mirrored at odd cell ordinals
   float fa = floor(a), ca = fa+1, fb = floor(b);
   vector lims = mirrorlimits((int)fa,a-fa,min(1,b-fa));
   // first possibly partial cell
   sum += ramp_sum_box(nl,0,lims.x,lims.y,bases,keys,vals);
   if(fb>fa) {
      // last possibly partial cell
      lims = mirrorlimits((int)fb,0,b-fb);
      sum += ramp_sum_box(nl,0,lims.x,lims.y,bases,keys,vals);
      // plus all the complete cells in between (there are 'nr' of them)
      int nr = (int)(fb-ca);
      if(nr>0) sum += nr * ramp_sum_box(nl,0,0,1,bases,keys,vals);
   }
   return sum;
}

vector ramp_summirror_box( SUM_BOX_PARMS(vector) ) {
   vector sum = 0;
   float fa = floor(a), ca = fa+1, fb = floor(b);
   vector lims = mirrorlimits((int)fa,a-fa,min(1,b-fa));
   sum += ramp_sum_box(nl,VZERO,lims.x,lims.y,bases,keys,vals);
   if(fb>fa) {
      lims = mirrorlimits((int)fb,0,b-fb);
      sum += ramp_sum_box(nl,VZERO,lims.x,lims.y,bases,keys,vals);
      int nr = (int)(fb-ca);
      if(nr>0) sum += nr * ramp_sum_box(nl,VZERO,0,1,bases,keys,vals);
   }
   return sum;
}

// Accum
float  ramp_sumaccum_box( SUM_BOX_PARMS(float ) ) {
   float sum = 0;
   // Same as cycle, except with a cummulative offset 'dv' per cell
   float fa = floor(a), ca = fa+1, fb = floor(b);
   float x = a-fa, y = min(1,b-fa), len = y-x;
   // first possibly partial cell
   sum += dv*fa*len + ramp_sum_box(nl,0,x,y,bases,keys,vals);
   if(fb>fa) {
      // last possibly partial cell
      len = b-fb;
      sum += dv*fb*len + ramp_sum_box(nl,0,0,b-fb,bases,keys,vals);
      // and the cells in between (there are 'nr' of them).
      // => nr*(sum of entire ramp) + dv*(sum of ints, or cell ordinals,
      //    between first and last complete cell)
      int nr = (int)(fb-ca);
      if(nr>0) sum += nr* ( ramp_sum_box(nl,0,0,1,bases,keys,vals) +
                            dv*(ca+fb-1)*0.5 );
   }
   return sum;
}

vector ramp_sumaccum_box( SUM_BOX_PARMS(vector) ) {
   vector sum = 0;
   float fa = floor(a), ca = fa+1, fb = floor(b);
   float x = a-fa, y = min(1,b-fa), len = y-x;
   sum += dv*fa*len + ramp_sum_box(nl,VZERO,x,y,bases,keys,vals);
   if(fb>fa) {
      len = b-fb;
      sum += dv*fb*len + ramp_sum_box(nl,VZERO,0,b-fb,bases,keys,vals);
      int nr = (int)(fb-ca);
      if(nr>0) sum += nr* ( ramp_sum_box(nl,VZERO,0,1,bases,keys,vals) +
                            dv*(ca+fb-1)*0.5 );
   }
   return sum;
}


#define RAMP_BOX_HEAD(T) \
   if(b<=a) return 0; \
   T  sum  = 0

// Float
float ramp_hold_hold_box( RAMP_FILTERED_PARMS(float) ) {
   RAMP_BOX_HEAD(float);
   if(a<0)        sum += ramp_sumconst_box(lv,min(b,0)-a);
   if(a<1 && b>0) sum += ramp_sum_box(nl,0,max(a,0),min(b,1),bases,keys,vals);
   if(b>1)        sum += ramp_sumconst_box(rv,b-max(a,1));
   return sum/(b-a);
}
float ramp_hold_slope_box( RAMP_FILTERED_PARMS(float) ) {
   RAMP_BOX_HEAD(float);
   if(a<0)        sum += ramp_sumconst_box(lv,min(b,0)-a);
   if(a<1 && b>0) sum += ramp_sum_box(nl,0,max(a,0),min(b,1),bases,keys,vals);
   if(b>1)        sum += ramp_sumslope_box(rv,rm,max(a,1)-1,b-1);
   return sum/(b-a);
}
float ramp_hold_cycle_box( RAMP_FILTERED_PARMS(float) ) {
   RAMP_BOX_HEAD(float);
   if(a<0) sum += ramp_sumconst_box(lv,min(b,0)-a);
   if(b>0) sum += ramp_sumcycle_box( SUM_BOX_ARGS(float,max(0,a),b) );
   return sum/(b-a);
}
float ramp_hold_mirror_box( RAMP_FILTERED_PARMS(float) ) {
   RAMP_BOX_HEAD(float);
   if(a<0) sum += ramp_sumconst_box(lv,min(b,0)-a);
   if(b>0) sum += ramp_summirror_box( SUM_BOX_ARGS(float,max(0,a),b) );
   return sum/(b-a);
}
float ramp_hold_accum_box( RAMP_FILTERED_PARMS(float) ) {
   RAMP_BOX_HEAD(float);
   if(a<0) sum += ramp_sumconst_box(lv,min(b,0)-a);
   if(b>0) sum += ramp_sumaccum_box( SUM_BOX_ARGS(float,max(0,a),b) );
   return sum/(b-a);
}

float ramp_slope_hold_box( RAMP_FILTERED_PARMS(float) ) {
   RAMP_BOX_HEAD(float);
   if(a<0)        sum += ramp_sumslope_box(lv,lm,a,min(0,b));
   if(a<1 && b>0) sum += ramp_sum_box(nl,0,max(a,0),min(b,1),bases,keys,vals);
   if(b>1)        sum += ramp_sumconst_box(rv,b-max(a,1));
   return sum/(b-a);
}
float ramp_slope_slope_box( RAMP_FILTERED_PARMS(float) ) {
   RAMP_BOX_HEAD(float);
   if(a<0)        sum += ramp_sumslope_box(lv,lm,a,min(0,b));
   if(a<1 && b>0) sum += ramp_sum_box(nl,0,max(a,0),min(b,1),bases,keys,vals);
   if(b>1)        sum += ramp_sumslope_box(rv,rm,max(a,1)-1,b-1);
   return sum/(b-a);
}
float ramp_slope_cycle_box( RAMP_FILTERED_PARMS(float) ) {
   RAMP_BOX_HEAD(float);
   if(a<0) sum += ramp_sumslope_box(lv,lm,a,min(0,b));
   if(b>0) sum += ramp_sumcycle_box( SUM_BOX_ARGS(float,max(0,a),b) );
   return sum/(b-a);
}
float ramp_slope_mirror_box( RAMP_FILTERED_PARMS(float) ) {
   RAMP_BOX_HEAD(float);
   if(a<0) sum += ramp_sumslope_box(lv,lm,a,min(0,b));
   if(b>0) sum += ramp_summirror_box( SUM_BOX_ARGS(float,max(0,a),b) );
   return sum/(b-a);
}
float ramp_slope_accum_box( RAMP_FILTERED_PARMS(float) ) {
   RAMP_BOX_HEAD(float);
   if(a<0) sum += ramp_sumslope_box(lv,lm,a,min(0,b));
   if(b>0) sum += ramp_sumaccum_box( SUM_BOX_ARGS(float,max(0,a),b) );
   return sum/(b-a);
}

float ramp_cycle_hold_box( RAMP_FILTERED_PARMS(float) ) {
   RAMP_BOX_HEAD(float);
   if(a<1) sum += ramp_sumcycle_box( SUM_BOX_ARGS(float,a,min(1,b)) );
   if(b>1) sum += ramp_sumconst_box(rv,b-max(a,1));
   return sum/(b-a);
}
float ramp_cycle_slope_box( RAMP_FILTERED_PARMS(float) ) {
   RAMP_BOX_HEAD(float);
   if(a<1) sum += ramp_sumcycle_box( SUM_BOX_ARGS(float,a,min(1,b)) );
   if(b>1) sum += ramp_sumslope_box(rv,rm,max(a,1)-1,b-1);
   return sum/(b-a);
}
float ramp_cycle_cycle_box( RAMP_FILTERED_PARMS(float) ) {
   RAMP_BOX_HEAD(float);
   sum += ramp_sumcycle_box( SUM_BOX_ARGS(float,a,b) );
   return sum/(b-a);
}
float ramp_cycle_mirror_box( RAMP_FILTERED_PARMS(float) ) {
   RAMP_BOX_HEAD(float);
   if(a<1) sum += ramp_sumcycle_box( SUM_BOX_ARGS(float,a,min(1,b)) );
   if(b>1) sum += ramp_summirror_box( SUM_BOX_ARGS(float,max(1,a),b) );
   return sum/(b-a);
}
float ramp_cycle_accum_box( RAMP_FILTERED_PARMS(float) ) {
   RAMP_BOX_HEAD(float);
   if(a<1) sum += ramp_sumcycle_box( SUM_BOX_ARGS(float,a,min(1,b)) );
   if(b>1) sum += ramp_sumaccum_box( SUM_BOX_ARGS(float,max(1,a),b) );
   return sum/(b-a);
}

float ramp_mirror_hold_box( RAMP_FILTERED_PARMS(float) ) {
   RAMP_BOX_HEAD(float);
   if(a<1) sum += ramp_summirror_box( SUM_BOX_ARGS(float,a,min(1,b)) );
   if(b>1) sum += ramp_sumconst_box(rv,b-max(a,1));
   return sum/(b-a);
}
float ramp_mirror_slope_box( RAMP_FILTERED_PARMS(float) ) {
   RAMP_BOX_HEAD(float);
   if(a<1) sum += ramp_summirror_box( SUM_BOX_ARGS(float,a,min(1,b)) );
   if(b>1) sum += ramp_sumslope_box(rv,rm,max(a,1)-1,b-1);
   return sum/(b-a);
}
float ramp_mirror_cycle_box( RAMP_FILTERED_PARMS(float) ) {
   RAMP_BOX_HEAD(float);
   if(a<1) sum += ramp_summirror_box( SUM_BOX_ARGS(float,a,min(1,b)) );
   if(b>1) sum += ramp_sumcycle_box( SUM_BOX_ARGS(float,max(1,a),b) );
   return sum/(b-a);
}
float ramp_mirror_mirror_box( RAMP_FILTERED_PARMS(float) ) {
   RAMP_BOX_HEAD(float);
   sum += ramp_summirror_box( SUM_BOX_ARGS(float,a,b) );
   return sum/(b-a);
}
float ramp_mirror_accum_box( RAMP_FILTERED_PARMS(float) ) {
   RAMP_BOX_HEAD(float);
   if(a<1) sum += ramp_summirror_box( SUM_BOX_ARGS(float,a,min(1,b)) );
   if(b>1) sum += ramp_sumaccum_box( SUM_BOX_ARGS(float,max(1,a),b) );
   return sum/(b-a);
}

float ramp_accum_hold_box( RAMP_FILTERED_PARMS(float) ) {
   RAMP_BOX_HEAD(float);
   if(a<1) sum += ramp_sumaccum_box( SUM_BOX_ARGS(float,a,min(1,b)) );
   if(b>1) sum += ramp_sumconst_box(rv,b-max(a,1));
   return sum/(b-a);
}
float ramp_accum_slope_box( RAMP_FILTERED_PARMS(float) ) {
   RAMP_BOX_HEAD(float);
   if(a<1) sum += ramp_sumaccum_box( SUM_BOX_ARGS(float,a,min(1,b)) );
   if(b>1) sum += ramp_sumslope_box(rv,rm,max(a,1)-1,b-1);
   return sum/(b-a);
}
float ramp_accum_cycle_box( RAMP_FILTERED_PARMS(float) ) {
   RAMP_BOX_HEAD(float);
   if(a<1) sum += ramp_sumaccum_box( SUM_BOX_ARGS(float,a,min(1,b)) );
   if(b>1) sum += ramp_sumcycle_box( SUM_BOX_ARGS(float,max(1,a),b) );
   return sum/(b-a);
}
float ramp_accum_mirror_box( RAMP_FILTERED_PARMS(float) ) {
   RAMP_BOX_HEAD(float);
   if(a<1) sum += ramp_sumaccum_box( SUM_BOX_ARGS(float,a,min(1,b)) );
   if(b>1) sum += ramp_summirror_box( SUM_BOX_ARGS(float,max(1,a),b) );
   return sum/(b-a);
}
float ramp_accum_accum_box( RAMP_FILTERED_PARMS(float) ) {
   RAMP_BOX_HEAD(float);
   sum += ramp_sumaccum_box( SUM_BOX_ARGS(float,a,b) );
   return sum/(b-a);
}


// Vector
vector ramp_hold_hold_box( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_BOX_HEAD(vector);
   if(a<0)        sum += ramp_sumconst_box(lv,min(b,0)-a);
   if(a<1 && b>0) sum += ramp_sum_box(nl,VZERO,max(a,0),min(b,1),bases,keys,vals);
   if(b>1)        sum += ramp_sumconst_box(rv,b-max(a,1));
   return sum/(b-a);
}
vector ramp_hold_slope_box( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_BOX_HEAD(vector);
   if(a<0)        sum += ramp_sumconst_box(lv,min(b,0)-a);
   if(a<1 && b>0) sum += ramp_sum_box(nl,VZERO,max(a,0),min(b,1),bases,keys,vals);
   if(b>1)        sum += ramp_sumslope_box(rv,rm,max(a,1)-1,b-1);
   return sum/(b-a);
}
vector ramp_hold_cycle_box( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_BOX_HEAD(vector);
   if(a<0) sum += ramp_sumconst_box(lv,min(b,0)-a);
   if(b>0) sum += ramp_sumcycle_box( SUM_BOX_ARGS(vector,max(0,a),b) );
   return sum/(b-a);
}
vector ramp_hold_mirror_box( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_BOX_HEAD(vector);
   if(a<0) sum += ramp_sumconst_box(lv,min(b,0)-a);
   if(b>0) sum += ramp_summirror_box( SUM_BOX_ARGS(vector,max(0,a),b) );
   return sum/(b-a);
}
vector ramp_hold_accum_box( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_BOX_HEAD(vector);
   if(a<0) sum += ramp_sumconst_box(lv,min(b,0)-a);
   if(b>0) sum += ramp_sumaccum_box( SUM_BOX_ARGS(vector,max(0,a),b) );
   return sum/(b-a);
}

vector ramp_slope_hold_box( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_BOX_HEAD(vector);
   if(a<0)        sum += ramp_sumslope_box(lv,lm,a,min(0,b));
   if(a<1 && b>0) sum += ramp_sum_box(nl,VZERO,max(a,0),min(b,1),bases,keys,vals);
   if(b>1)        sum += ramp_sumconst_box(rv,b-max(a,1));
   return sum/(b-a);
}
vector ramp_slope_slope_box( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_BOX_HEAD(vector);
   if(a<0)        sum += ramp_sumslope_box(lv,lm,a,min(0,b));
   if(a<1 && b>0) sum += ramp_sum_box(nl,VZERO,max(a,0),min(b,1),bases,keys,vals);
   if(b>1)        sum += ramp_sumslope_box(rv,rm,max(a,1)-1,b-1);
   return sum/(b-a);
}
vector ramp_slope_cycle_box( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_BOX_HEAD(vector);
   if(a<0) sum += ramp_sumslope_box(lv,lm,a,min(0,b));
   if(b>0) sum += ramp_sumcycle_box( SUM_BOX_ARGS(vector,max(0,a),b) );
   return sum/(b-a);
}
vector ramp_slope_mirror_box( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_BOX_HEAD(vector);
   if(a<0) sum += ramp_sumslope_box(lv,lm,a,min(0,b));
   if(b>0) sum += ramp_summirror_box( SUM_BOX_ARGS(vector,max(0,a),b) );
   return sum/(b-a);
}
vector ramp_slope_accum_box( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_BOX_HEAD(vector);
   if(a<0) sum += ramp_sumslope_box(lv,lm,a,min(0,b));
   if(b>0) sum += ramp_sumaccum_box( SUM_BOX_ARGS(vector,max(0,a),b) );
   return sum/(b-a);
}

vector ramp_cycle_hold_box( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_BOX_HEAD(vector);
   if(a<1) sum += ramp_sumcycle_box( SUM_BOX_ARGS(vector,a,min(1,b)) );
   if(b>1) sum += ramp_sumconst_box(rv,b-max(a,1));
   return sum/(b-a);
}
vector ramp_cycle_slope_box( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_BOX_HEAD(vector);
   if(a<1) sum += ramp_sumcycle_box( SUM_BOX_ARGS(vector,a,min(1,b)) );
   if(b>1) sum += ramp_sumslope_box(rv,rm,max(a,1)-1,b-1);
   return sum/(b-a);
}
vector ramp_cycle_cycle_box( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_BOX_HEAD(vector);
   sum += ramp_sumcycle_box( SUM_BOX_ARGS(vector,a,b) );
   return sum/(b-a);
}
vector ramp_cycle_mirror_box( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_BOX_HEAD(vector);
   if(a<1) sum += ramp_sumcycle_box( SUM_BOX_ARGS(vector,a,min(1,b)) );
   if(b>1) sum += ramp_summirror_box( SUM_BOX_ARGS(vector,max(1,a),b) );
   return sum/(b-a);
}
vector ramp_cycle_accum_box( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_BOX_HEAD(vector);
   if(a<1) sum += ramp_sumcycle_box( SUM_BOX_ARGS(vector,a,min(1,b)) );
   if(b>1) sum += ramp_sumaccum_box( SUM_BOX_ARGS(vector,max(1,a),b) );
   return sum/(b-a);
}

vector ramp_mirror_hold_box( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_BOX_HEAD(vector);
   if(a<1) sum += ramp_summirror_box( SUM_BOX_ARGS(vector,a,min(1,b)) );
   if(b>1) sum += ramp_sumconst_box(rv,b-max(a,1));
   return sum/(b-a);
}
vector ramp_mirror_slope_box( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_BOX_HEAD(vector);
   if(a<1) sum += ramp_summirror_box( SUM_BOX_ARGS(vector,a,min(1,b)) );
   if(b>1) sum += ramp_sumslope_box(rv,rm,max(a,1)-1,b-1);
   return sum/(b-a);
}
vector ramp_mirror_cycle_box( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_BOX_HEAD(vector);
   if(a<1) sum += ramp_summirror_box( SUM_BOX_ARGS(vector,a,min(1,b)) );
   if(b>1) sum += ramp_sumcycle_box( SUM_BOX_ARGS(vector,max(1,a),b) );
   return sum/(b-a);
}
vector ramp_mirror_mirror_box( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_BOX_HEAD(vector);
   sum += ramp_summirror_box( SUM_BOX_ARGS(vector,a,b) );
   return sum/(b-a);
}
vector ramp_mirror_accum_box( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_BOX_HEAD(vector);
   if(a<1) sum += ramp_summirror_box( SUM_BOX_ARGS(vector,a,min(1,b)) );
   if(b>1) sum += ramp_sumaccum_box( SUM_BOX_ARGS(vector,max(1,a),b) );
   return sum/(b-a);
}

vector ramp_accum_hold_box( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_BOX_HEAD(vector);
   if(a<1) sum += ramp_sumaccum_box( SUM_BOX_ARGS(vector,a,min(1,b)) );
   if(b>1) sum += ramp_sumconst_box(rv,b-max(a,1));
   return sum/(b-a);
}
vector ramp_accum_slope_box( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_BOX_HEAD(vector);
   if(a<1) sum += ramp_sumaccum_box( SUM_BOX_ARGS(vector,a,min(1,b)) );
   if(b>1) sum += ramp_sumslope_box(rv,rm,max(a,1)-1,b-1);
   return sum/(b-a);
}
vector ramp_accum_cycle_box( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_BOX_HEAD(vector);
   if(a<1) sum += ramp_sumaccum_box( SUM_BOX_ARGS(vector,a,min(1,b)) );
   if(b>1) sum += ramp_sumcycle_box( SUM_BOX_ARGS(vector,max(1,a),b) );
   return sum/(b-a);
}
vector ramp_accum_mirror_box( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_BOX_HEAD(vector);
   if(a<1) sum += ramp_sumaccum_box( SUM_BOX_ARGS(vector,a,min(1,b)) );
   if(b>1) sum += ramp_summirror_box( SUM_BOX_ARGS(vector,max(1,a),b) );
   return sum/(b-a);
}
vector ramp_accum_accum_box( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_BOX_HEAD(vector);
   sum += ramp_sumaccum_box( SUM_BOX_ARGS(vector,a,b) );
   return sum/(b-a);
}






// FILTER: Gaussian
//--------------------------------------
#define RAMP_GAUSS_NORM  0.313308687321307169163532
#define RAMP_GAUSS_INORM 3.191740415976627300594819

// native width of this filter
float fscale_gauss() { return 1.25; }

// convolve with a constant -- for "hold" boundaries
float ramp_sumconst_gauss(float val,fa,fb) {
   float erfa = erf(4*C_SQRT2*fa);
   float erfb = erf(4*C_SQRT2*fb);
   return 0.125 * val * C_SQRTPI_2 * (erfb-erfa);
}
vector ramp_sumconst_gauss(vector val; float fa,fb) {
   float erfa = erf(4*C_SQRT2*fa);
   float erfb = erf(4*C_SQRT2*fb);
   return 0.125 * val * C_SQRTPI_2 * (erfb-erfa);
}

// convolve with a line (mx+b) -- for "slope" boundaries
// computed for key at position 'pos' with value 'val' and slope 'slope'
// convolved with gaussian centered at 'c' of width 'w'
float ramp_sumslope_gauss(float pos,val,slope,fa,fb, c,w) {
   if(slope==0) return ramp_sumconst_gauss(val,fa,fb);
   float b    = -slope*pos + val;
   float def  = exp(-32*fa*fa)-exp(-32*fb*fb);
   float derf = erf(4*C_SQRT2*fb) - erf(4*C_SQRT2*fa);
   return (1/64)*(def*slope*w+4*derf*(b+c*slope)*C_SQRT2PI);
}
vector ramp_sumslope_gauss(float pos; vector val,slope; float fa,fb,c,w) {
   if(slope==VZERO) return ramp_sumconst_gauss(val,fa,fb);
   vector b    = -slope*pos + val;
   float  def  = exp(-32*fa*fa)-exp(-32*fb*fb);
   float  derf = erf(4*C_SQRT2*fb) - erf(4*C_SQRT2*fa);
   return (1/64)*(def*slope*w+4*derf*(b+c*slope)*C_SQRT2PI);
}


// Convolve gaussian with a single ramp segment
float ramp_sumseg_gauss(
      int    n;         // array index at segment start
      int    nl;        // last array index
      float  ci;        // const of integration
      float  sa,sb;     // in-segment limits of integration ([0,1], a<=b)
      float  fa,fb;     // in-filter limits of integration ([-0.5,+0.5], fa<=fb)
      float  c,w,iw;    // center, scale, and inverse-scale of gaussian filter
      string bases[];   // ramp basis strings
      float  keys[];    // ramp key positions
      float  vals[]     // ramp key values
   ) 
{
   // Convolution done in filter space (transformed segment) as opposed to
   // segment space (transformed filter) for reasons of numerical stability.
   // In filter space, the exponential arguments are guaranteed to be in
   // (+/-)32*(fa^2+fb^2), or (+/-32)*[0,0.5] -> an absolute max of 16, and
   // exp(+/-16) is still safe. In segment space this would blow up immediately
   // due to the unbounded exponentials in the varying-width gaussian.
   // All calculations assume: w>0 && fa<fb in [-0.5,0.5] && sa<sb in [0,1]

   float erfa = erf(4*C_SQRT2*fa);
   float erfb = erf(4*C_SQRT2*fb);
   float derf = erfb-erfa;

   float dfab  = fa - fb;
   float adfab = abs(dfab);

   // If the size of the filter-space integration segment 'abs(fb-fa)'is below a 
   // threshold where differencing will start to underflow (around 8.7e-5, which 
   // is equivalent to a filter size of ~6,000 on a unit segment), then we just 
   // return a single weighted sample at the center. Note: this issue only shows
   // when the center of the kernel is at or close to w units away from either
   // end of the segment and only when w is very large
   if(adfab<8.7e-5)
      return 0.125 * (0.5*(vals[n]+vals[n+1])+ci) * C_SQRTPI_2 * derf;

   float len = keys[n+1]-keys[n];
   int   n0  = max(0,n-1), n1 = n, n2 = min(nl,n+1), n3 = min(nl,n2+1);
   float v0  = vals[n0]+ci, v1 = vals[n1]+ci, 
         v2  = vals[n2]+ci, v3 = vals[n3]+ci;

   float fa2 = fa*fa, fa3 = fa2*fa;
   float fb2 = fb*fb, fb3 = fb2*fb;
   float efa = exp(32*fa2), nefa = 1/efa;
   float efb = exp(32*fb2), nefb = 1/efb;
   float efab = efa*efb,    nefab = nefa*nefb;

   if(bases[n]=="constant") {
      return 0.125 * v1 * C_SQRTPI_2 * derf;
   }

   // If the size of the filter-space integration segment is below a threshold 
   // where the approximation to the erf() function starts to become unstable 
   // (around 0.006 for the implementation above), then all interpolants are 
   // treated as linear (with only a loss of precission in the 1e-5 region for 
   // the cubics, which is acceptable for this application).

   float m,b;

   if(bases[n]=="linear" || adfab<0.055) {      // mx + b
      m  = (v2-v1) / len;
      b  = -m*keys[n] + v1;
      return ((nefa-nefb)*m + 4*iw*C_SQRT2PI*(c*m + b)*derf) / (iw*64.0);
   }

   // The transform of x [fa,fb] from filter space to segment space is:
   // sa - (fa*(sb-sa) + x*(sb-sa)) / (fb-fa)
   // but the convolution is calculated using the simplest form (m*x + b) to
   // reduce the otherwise hundreds of terms. So in that form, the scale 'm' 
   // and offset 'b' become:  m = (sb-sa) / (fb-fa), and b = sa - m*fa
   float m2,m3,b2,b3;
   m   = (sb-sa)/(fb-fa); m2 = m*m;   m3 = m2*m;
   b   = sa - m*fa;       b2 = b*b;   b3 = b2*b;
   float bm1 = b-1,       bm12 = bm1*bm1;

   if(bases[n]=="catmull-rom")
   {
      // TODO: collect w.r.t: v0,v1,v2,v3, then each term w.r.t: m,b, 
      // then back to matrix form if possible
      float iderf = erfa-erfb;
      return 
         (nefab/4096)*(
            (-efb*(1+32*fa2)+efa*(1+32*fb2))*m3*(v0-3*v1+3*v2-v3)
           +32*(efa-efb)*m*(
              bm1*(-1+3*b)*v0-v2+2*b*(5*v1-4*v2+v3)-3*b2*(3*v1-3*v2+v3))
           +2*m2*((-2+3*b)*v0+5*v1-4*v2+v3-3*b*(3*v1-3*v2+v3))*(
              16*(-efb*fa+efa*fb)+efab*C_SQRT2PI*iderf)
           +128*efab*C_SQRT2PI*iderf*(
              -2*v1+b*(bm12*v0-v2+b*((5-3*b)*v1-4*v2+3*b*v2+v3-b*v3)))
         );
   }

   else if(bases[n]=="monotonecubic")
   {
      float m0 = ramp_slope_mc(v0,v1,v2);
      float m1 = ramp_slope_mc(v1,v2,v3);

      // TODO: collect w.r.t: v1,v2,m0,m1, then each term w.r.t: m,b, 
      // then back to matrix form if possible
      return
         (nefab/2048)*(
            -32*(efa-efb)*m*(bm1*(-1+3*b)*m0+b*((-2+3*b)*m1+6*bm1*(v1-v2)))
            -(-efb*(1+32*fa2)+efa*(1+32*fb2))*m3*(m0+m1+2*v1-2*v2)
            -2*m2*((-2+3*b)*m0+(-1+3*b)*m1+3*(-1+2*b)*(v1-v2))*(
               16*(-efb*fa+efa*fb)+efab*C_SQRT2PI*(erfa-erfb))
            +128*efab*C_SQRT2PI*derf*(v1+b*(bm12*m0+b*(bm1*m1+(-3+2*b)*(v1-v2))))
         );

   }

   return 0;

}
vector ramp_sumseg_gauss(
      int    n;         // array index at segment start
      int    nl;        // last array index
      vector ci;        // const of integration
      float  sa,sb;     // in-segment limits of integration ([0,1], a<=b)
      float  fa,fb;     // in-filter limits of integration ([-0.5,+0.5], fa<=fb)
      float  c,w,iw;    // center, scale, and inverse-scale of gaussian filter
      string bases[];   // ramp basis strings
      float  keys[];    // ramp key positions
      vector vals[]     // ramp key values
   ) 
{
   // See the float version above for comments
   float erfa = erf(4*C_SQRT2*fa);
   float erfb = erf(4*C_SQRT2*fb);
   float derf = erfb-erfa;

   float dfab  = fa - fb;
   float adfab = abs(dfab);

   if(adfab<8.7e-5)
      return 0.125 * 0.5*(vals[n]+vals[n+1]) * C_SQRTPI_2 * derf;

   float  len = keys[n+1]-keys[n];
   int    n0  = max(0,n-1), n1 = n, n2 = min(nl,n+1), n3 = min(nl,n2+1);
   vector v0  = vals[n0]+ci, v1 = vals[n1]+ci, 
          v2  = vals[n2]+ci, v3 = vals[n3]+ci;

   float fa2 = fa*fa, fa3 = fa2*fa;
   float fb2 = fb*fb, fb3 = fb2*fb;
   float efa = exp(32*fa2), nefa = 1/efa;
   float efb = exp(32*fb2), nefb = 1/efb;
   float efab = efa*efb,    nefab = nefa*nefb;

   if(bases[n]=="constant") {
      return 0.125 * v1 * C_SQRTPI_2 * derf;
   }

   vector m,b;

   if(bases[n]=="linear" || adfab<0.015) {      // mx + b
      m  = (v2-v1) / len;
      b  = -m*keys[n] + v1;
      return ((nefa-nefb)*m + 4*iw*C_SQRT2PI*(c*m + b)*derf) / (iw*64.0);
   }

   vector m2,m3,b2,b3;
   m   = (sb-sa)/(fb-fa);   m2 = m*m;  m3 = m2*m;
   b   = (vector)sa - m*fa; b2 = b*b;  b3 = b2*b;
   vector bm1 = b-1, bm12 = bm1*bm1;

   if(bases[n]=="catmull-rom")
   {
      float iderf = erfa-erfb;
      return 
         (nefab/4096)*(
            (-efb*(1+32*fa2)+efa*(1+32*fb2))*m3*(v0-3*v1+3*v2-v3)
           +32*(efa-efb)*m*(
              bm1*(-1+3*b)*v0-v2+2*b*(5*v1-4*v2+v3)-3*b2*(3*v1-3*v2+v3))
           +2*m2*((-2+3*b)*v0+5*v1-4*v2+v3-3*b*(3*v1-3*v2+v3))*(
              16*(-efb*fa+efa*fb)+efab*C_SQRT2PI*iderf)
           +128*efab*C_SQRT2PI*iderf*(
              -2*v1+ b*(bm12*v0-v2+b*(((vector)5-3*b)*v1-4*v2+3*b*v2+v3-b*v3)))
         );
   }

   else if(bases[n]=="monotonecubic")
   {
      vector m0 = ramp_slope_mc(v0,v1,v2);
      vector m1 = ramp_slope_mc(v1,v2,v3);

      return
         (nefab/2048)*(
            -32*(efa-efb)*m*(bm1*(-1+3*b)*m0+b*((-2+3*b)*m1+6*bm1*(v1-v2)))
            -(-efb*(1+32*fa2)+efa*(1+32*fb2))*m3*(m0+m1+2*v1-2*v2)
            -2*m2*((-2+3*b)*m0+(-1+3*b)*m1+3*(-1+2*b)*(v1-v2))*(
               16*(-efb*fa+efa*fb)+efab*C_SQRT2PI*(erfa-erfb))
            +128*efab*C_SQRT2PI*derf*(v1+b*(bm12*m0+b*(bm1*m1+(-3+2*b)*(v1-v2))))
         );

   }

   return 0;

}

// convolve gaussian with a portion {a,b} of the entire [0,1] ramp
float ramp_sum_gauss (
      int    nl;        // last array index
      float  ci;        // const of integration
      float  a,b;       // real limits of integration in [0,1]
      float  c,w,iw;    // center, scale, and inverse-scale of gaussian filter
      string bases[];   // ramp basis strings
      float  keys[];    // ramp key positions
      float  vals[];    // ramp key values
   ) 
{
   float sum = 0;

   float f0 = (spline("solvelinear",a,keys));
   float f1 = (spline("solvelinear",b,keys));
   int   n0 = floor(f0*nl);
   int   n1 = (int)ceil(f1*nl);
   int n; float k0,k1,ka,kb,il,sa,sb,fa,fb;
   for(n=n0;n<n1;n++) { 
      k0 = keys[n];
      k1 = keys[n+1];
      il = 1/(k1-k0);
      ka = max(k0,a);
      kb = min(k1,b);
      sa = (ka-k0)*il;  // segment limits {sa,sb} in [0,1]
      sb = (kb-k0)*il;
      fa = (ka-c) *iw;  // filter limits {fa,fb} in [-0.5,0.5]
      fb = (kb-c) *iw;
      sum += ramp_sumseg_gauss(n,nl, ci, sa,sb, fa,fb, c,w,iw, bases,keys,vals);
   }

   return sum;
}
vector ramp_sum_gauss (
      int    nl;        // last array index
      vector ci;        // const of integration
      float  a,b;       // limits of integration in [0,1]
      float  c,w,iw;    // center, scale, and inverse-scale of gaussian filter
      string bases[];   // ramp basis strings
      float  keys[];    // ramp key positions
      vector vals[];    // ramp key values
   ) 
{
   vector sum = 0;

   float f0 = (spline("solvelinear",a,keys));
   float f1 = (spline("solvelinear",b,keys));
   int   n0 = (int)float(floor(f0*nl));
   int   n1 = (int)float(ceil(f1*nl));
   int n; float k0,k1,ka,kb,il,sa,sb,fa,fb;
   for(n=n0;n<n1;n++) { 
      k0 = keys[n];
      k1 = keys[n+1];
      il = 1/(k1-k0);
      ka = max(k0,a);
      kb = min(k1,b);
      sa = (ka-k0)*il;  // segment limits {sa,sb} in [0,1]
      sb = (kb-k0)*il;
      fa = (ka-c) *iw;  // filter limits {fa,fb} in [-0.5,0.5]
      fb = (kb-c) *iw;
      sum += ramp_sumseg_gauss(n,nl, ci, sa,sb, fa,fb, c,w,iw, bases,keys,vals);
   }

   return sum;
}


#define SUM_GAUSS_PARMS(T) \
   int nl; float a,b,x,w,iw; T dv; string bases[]; float keys[]; T vals[];
#define SUM_GAUSS_ARGS(T,a,b,x,w) \
   nl, a,b,x,w,iw, dv, bases, keys, vals

// Cycle
float  ramp_sumcycle_gauss(SUM_GAUSS_PARMS(float )) {
   float sum = 0;
   float xk  = floor(a);
   vector xv;
   while(xk<b) {
      xv = set(max(xk,a),min(xk+1,b),x) - xk;
      sum += ramp_sum_gauss(nl,0,xv[0],xv[1],xv[2],w,iw,bases,keys,vals);
      xk += 1;
   }
   return sum;
}
vector ramp_sumcycle_gauss( SUM_GAUSS_PARMS(vector) ) {
   vector sum = 0;
   float xk  = floor(a);
   vector xv;
   while(xk<b) {
      xv = set(max(xk,a),min(xk+1,b),x) - xk;
      sum += ramp_sum_gauss(nl,VZERO,xv[0],xv[1],xv[2],w,iw,bases,keys,vals);
      xk += 1;
   }
   return sum;
}

// Mirror
float  ramp_summirror_gauss( SUM_GAUSS_PARMS(float ) ) {
   float  sum = 0;
   float  xk  = floor(a);
   vector xv[]; resize(xv,2);
   int nm;
   while(xk<b) {
      nm = ((int)xk)%2;
      xv[0] = set(max(xk,a),min(xk+1,b),x) - xk;
      xv[1] = VONE - set(xv[0].y,xv[0].x,xv[0].z);
      sum += ramp_sum_gauss(nl,0,xv[nm].x,xv[nm].y,xv[nm].z,w,iw,bases,keys,vals);
      xk += 1;
   }
   return sum;
}
vector ramp_summirror_gauss( SUM_GAUSS_PARMS(vector) ) {
   vector sum = 0;
   float  xk  = floor(a);
   vector xv[]; resize(xv,2);
   int nm;
   while(xk<b) {
      nm = ((int)xk)%2;
      xv[0] = set(max(xk,a),min(xk+1,b),x) - xk;
      xv[1] = VONE - set(xv[0].y,xv[0].x,xv[0].z);
      sum += ramp_sum_gauss(nl,VZERO,xv[nm].x,xv[nm].y,xv[nm].z,w,iw,bases,keys,vals);
      xk += 1;
   }
   return sum;
}

// Accum
float  ramp_sumaccum_gauss( SUM_GAUSS_PARMS(float ) ) {
   float sum = 0;
   float xk  = floor(a - (a>0));
   vector xv;
   while(xk<b) {
      xv = set(max(xk,a),min(xk+1,b),x) - xk;
      sum += ramp_sum_gauss(nl,dv*xk,xv[0],xv[1],xv[2],w,iw,bases,keys,vals);
      xk += 1;
   }
   return sum;
}
vector ramp_sumaccum_gauss( SUM_GAUSS_PARMS(vector) ) {
   vector sum = 0;
   float  xk  = floor(a - (a>0));
   vector xv;
   while(xk<b) {
      xv = set(max(xk,a),min(xk+1,b),x) - xk;
      sum += ramp_sum_gauss(nl,dv*xk,xv[0],xv[1],xv[2],w,iw,bases,keys,vals);
      xk += 1;
   }
   return sum;
}


#define RAMP_GAUSS_HEAD(type) \
   if(b<=a) return 0; \
   float w    = (b-a); \
   float iw   = 1.0/w; \
   type  sum  = 0

// Float
float ramp_hold_hold_gauss( RAMP_FILTERED_PARMS(float) ) {
   RAMP_GAUSS_HEAD(float);
   if(a<0)        sum += ramp_sumconst_gauss(lv,-0.5,(min(b,0)-x)*iw);
   if(a<1 && b>0) sum += ramp_sum_gauss(nl,0,max(a,0),min(b,1),x,w,iw,bases,keys,vals);
   if(b>1)        sum += ramp_sumconst_gauss(rv,(max(a,1)-x)*iw,0.5);
   return sum*RAMP_GAUSS_INORM;
}
float ramp_hold_slope_gauss( RAMP_FILTERED_PARMS(float) ) {
   RAMP_GAUSS_HEAD(float);
   if(a<0)        sum += ramp_sumconst_gauss(lv,-0.5,(min(b,0)-x)*iw);
   if(a<1 && b>0) sum += ramp_sum_gauss(nl,0,max(a,0),min(b,1),x,w,iw,bases,keys,vals);
   if(b>1)        sum += ramp_sumslope_gauss(1,rv,rm,(max(a,1)-x)*iw,0.5,x,w);
   return sum*RAMP_GAUSS_INORM;
}
float ramp_hold_cycle_gauss( RAMP_FILTERED_PARMS(float) ) {
   RAMP_GAUSS_HEAD(float);
   if(a<0) sum += ramp_sumconst_gauss(lv,-0.5,(min(b,0)-x)*iw);
   if(b>0) sum += ramp_sumcycle_gauss(SUM_GAUSS_ARGS(float,max(0,a),b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}
float ramp_hold_mirror_gauss( RAMP_FILTERED_PARMS(float) ) {
   RAMP_GAUSS_HEAD(float);
   if(a<0) sum += ramp_sumconst_gauss(lv,-0.5,(min(b,0)-x)*iw);
   if(b>0) sum += ramp_summirror_gauss(SUM_GAUSS_ARGS(float,max(0,a),b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}
float ramp_hold_accum_gauss( RAMP_FILTERED_PARMS(float) ) {
   RAMP_GAUSS_HEAD(float);
   if(a<0) sum += ramp_sumconst_gauss(lv,-0.5,(min(b,0)-x)*iw);
   if(b>0) sum += ramp_sumaccum_gauss(SUM_GAUSS_ARGS(float,max(0,a),b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}

float ramp_slope_hold_gauss( RAMP_FILTERED_PARMS(float) ) {
   RAMP_GAUSS_HEAD(float);
   if(a<0)        sum += ramp_sumslope_gauss(0,lv,lm,-0.5,(min(b,0)-x)*iw,x,w);
   if(a<1 && b>0) sum += ramp_sum_gauss(nl,0,max(a,0),min(b,1),x,w,iw,bases,keys,vals);
   if(b>1)        sum += ramp_sumconst_gauss(rv,(max(a,1)-x)*iw,0.5);
   return sum*RAMP_GAUSS_INORM;
}
float ramp_slope_slope_gauss( RAMP_FILTERED_PARMS(float) ) {
   RAMP_GAUSS_HEAD(float);
   if(a<0)        sum += ramp_sumslope_gauss(0,lv,lm,-0.5,(min(b,0)-x)*iw,x,w);
   if(a<1 && b>0) sum += ramp_sum_gauss(nl,0,max(a,0),min(b,1),x,w,iw,bases,keys,vals);
   if(b>1)        sum += ramp_sumslope_gauss(1,rv,rm,(max(a,1)-x)*iw,0.5,x,w);
   return sum*RAMP_GAUSS_INORM;
}
float ramp_slope_cycle_gauss( RAMP_FILTERED_PARMS(float) ) {
   RAMP_GAUSS_HEAD(float);
   if(a<0)        sum += ramp_sumslope_gauss(0,lv,lm,-0.5,(min(b,0)-x)*iw,x,w);
   if(b>0) sum += ramp_sumcycle_gauss(SUM_GAUSS_ARGS(float,max(0,a),b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}
float ramp_slope_mirror_gauss( RAMP_FILTERED_PARMS(float) ) {
   RAMP_GAUSS_HEAD(float);
   if(a<0)        sum += ramp_sumslope_gauss(0,lv,lm,-0.5,(min(b,0)-x)*iw,x,w);
   if(b>0) sum += ramp_summirror_gauss(SUM_GAUSS_ARGS(float,max(0,a),b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}
float ramp_slope_accum_gauss( RAMP_FILTERED_PARMS(float) ) {
   RAMP_GAUSS_HEAD(float);
   if(a<0)        sum += ramp_sumslope_gauss(0,lv,lm,-0.5,(min(b,0)-x)*iw,x,w);
   if(b>0) sum += ramp_sumaccum_gauss(SUM_GAUSS_ARGS(float,max(0,a),b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}

float ramp_cycle_hold_gauss( RAMP_FILTERED_PARMS(float) ) {
   RAMP_GAUSS_HEAD(float);
   if(a<1) sum += ramp_sumcycle_gauss(SUM_GAUSS_ARGS(float,a,min(1,b),x,w) );
   if(b>1)        sum += ramp_sumconst_gauss(rv,(max(a,1)-x)*iw,0.5);
   return sum*RAMP_GAUSS_INORM;
}
float ramp_cycle_slope_gauss( RAMP_FILTERED_PARMS(float) ) {
   RAMP_GAUSS_HEAD(float);
   if(a<1) sum += ramp_sumcycle_gauss(SUM_GAUSS_ARGS(float,a,min(1,b),x,w) );
   if(b>1)        sum += ramp_sumslope_gauss(1,rv,rm,(max(a,1)-x)*iw,0.5,x,w);
   return sum*RAMP_GAUSS_INORM;
}
float ramp_cycle_cycle_gauss( RAMP_FILTERED_PARMS(float) ) {
   RAMP_GAUSS_HEAD(float);
   sum = ramp_sumcycle_gauss(SUM_GAUSS_ARGS(float,a,b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}
float ramp_cycle_mirror_gauss( RAMP_FILTERED_PARMS(float) ) {
   RAMP_GAUSS_HEAD(float);
   if(a<1) sum += ramp_sumcycle_gauss(SUM_GAUSS_ARGS(float,a,min(1,b),x,w) );
   if(b>1) sum += ramp_summirror_gauss(SUM_GAUSS_ARGS(float,max(1,a),b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}
float ramp_cycle_accum_gauss( RAMP_FILTERED_PARMS(float) ) {
   RAMP_GAUSS_HEAD(float);
   if(a<1) sum += ramp_sumcycle_gauss(SUM_GAUSS_ARGS(float,a,min(1,b),x,w) );
   if(b>1) sum += ramp_sumaccum_gauss(SUM_GAUSS_ARGS(float,max(1,a),b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}

float ramp_mirror_hold_gauss( RAMP_FILTERED_PARMS(float) ) {
   RAMP_GAUSS_HEAD(float);
   if(a<1) sum += ramp_summirror_gauss(SUM_GAUSS_ARGS(float,a,min(1,b),x,w) );
   if(b>1) sum += ramp_sumconst_gauss(rv,(max(a,1)-x)*iw,0.5);
   return sum*RAMP_GAUSS_INORM;
}
float ramp_mirror_slope_gauss( RAMP_FILTERED_PARMS(float) ) {
   RAMP_GAUSS_HEAD(float);
   if(a<1) sum += ramp_summirror_gauss(SUM_GAUSS_ARGS(float,a,min(1,b),x,w) );
   if(b>1) sum += ramp_sumslope_gauss(1,rv,rm,(max(a,1)-x)*iw,0.5,x,w);
   return sum*RAMP_GAUSS_INORM;
}
float ramp_mirror_cycle_gauss( RAMP_FILTERED_PARMS(float) ) {
   RAMP_GAUSS_HEAD(float);
   if(a<1) sum += ramp_summirror_gauss(SUM_GAUSS_ARGS(float,a,min(1,b),x,w) );
   if(b>1) sum += ramp_sumcycle_gauss(SUM_GAUSS_ARGS(float,max(1,a),b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}
float ramp_mirror_mirror_gauss( RAMP_FILTERED_PARMS(float) ) {
   RAMP_GAUSS_HEAD(float);
   sum += ramp_summirror_gauss(SUM_GAUSS_ARGS(float,a,b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}
float ramp_mirror_accum_gauss( RAMP_FILTERED_PARMS(float) ) {
   RAMP_GAUSS_HEAD(float);
   if(a<1) sum += ramp_summirror_gauss(SUM_GAUSS_ARGS(float,a,min(1,b),x,w) );
   if(b>1) sum += ramp_sumaccum_gauss(SUM_GAUSS_ARGS(float,max(1,a),b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}

float ramp_accum_hold_gauss( RAMP_FILTERED_PARMS(float) ) {
   RAMP_GAUSS_HEAD(float);
   if(a<1) sum += ramp_sumaccum_gauss(SUM_GAUSS_ARGS(float,a,min(1,b),x,w) );
   if(b>1) sum += ramp_sumconst_gauss(rv,(max(a,1)-x)*iw,0.5);
   return sum*RAMP_GAUSS_INORM;
}
float ramp_accum_slope_gauss( RAMP_FILTERED_PARMS(float) ) {
   RAMP_GAUSS_HEAD(float);
   if(a<1) sum += ramp_sumaccum_gauss(SUM_GAUSS_ARGS(float,a,min(1,b),x,w) );
   if(b>1) sum += ramp_sumslope_gauss(1,rv,rm,(max(a,1)-x)*iw,0.5,x,w);
   return sum*RAMP_GAUSS_INORM;
}
float ramp_accum_cycle_gauss( RAMP_FILTERED_PARMS(float) ) {
   RAMP_GAUSS_HEAD(float);
   if(a<1) sum += ramp_sumaccum_gauss(SUM_GAUSS_ARGS(float,a,min(1,b),x,w) );
   if(b>1) sum += ramp_sumcycle_gauss(SUM_GAUSS_ARGS(float,max(1,a),b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}
float ramp_accum_mirror_gauss( RAMP_FILTERED_PARMS(float) ) {
   RAMP_GAUSS_HEAD(float);
   if(a<1) sum += ramp_sumaccum_gauss(SUM_GAUSS_ARGS(float,a,min(1,b),x,w) );
   if(b>1) sum += ramp_summirror_gauss(SUM_GAUSS_ARGS(float,max(1,a),b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}
float ramp_accum_accum_gauss( RAMP_FILTERED_PARMS(float) ) {
   RAMP_GAUSS_HEAD(float);
   sum += ramp_sumaccum_gauss(SUM_GAUSS_ARGS(float,a,b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}



// Vector
vector ramp_hold_hold_gauss( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_GAUSS_HEAD(vector);
   if(a<0)        sum += ramp_sumconst_gauss(lv,-0.5,(min(b,0)-x)*iw);
   if(a<1 && b>0) 
      sum += ramp_sum_gauss(nl,VZERO,max(a,0),min(b,1),x,w,iw,bases,keys,vals);
   if(b>1)        sum += ramp_sumconst_gauss(rv,(max(a,1)-x)*iw,0.5);
   return sum*RAMP_GAUSS_INORM;
}
vector ramp_hold_slope_gauss( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_GAUSS_HEAD(vector);
   if(a<0)        sum += ramp_sumconst_gauss(lv,-0.5,(min(b,0)-x)*iw);
   if(a<1 && b>0) 
      sum += ramp_sum_gauss(nl,VZERO,max(a,0),min(b,1),x,w,iw,bases,keys,vals);
   if(b>1)        sum += ramp_sumslope_gauss(1,rv,rm,(max(a,1)-x)*iw,0.5,x,w);
   return sum*RAMP_GAUSS_INORM;
}
vector ramp_hold_cycle_gauss( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_GAUSS_HEAD(vector);
   if(a<0) sum += ramp_sumconst_gauss(lv,-0.5,(min(b,0)-x)*iw);
   if(b>0) sum += ramp_sumcycle_gauss(SUM_GAUSS_ARGS(vector,max(0,a),b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}
vector ramp_hold_mirror_gauss( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_GAUSS_HEAD(vector);
   if(a<0) sum += ramp_sumconst_gauss(lv,-0.5,(min(b,0)-x)*iw);
   if(b>0) sum += ramp_summirror_gauss(SUM_GAUSS_ARGS(vector,max(0,a),b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}
vector ramp_hold_accum_gauss( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_GAUSS_HEAD(vector);
   if(a<0) sum += ramp_sumconst_gauss(lv,-0.5,(min(b,0)-x)*iw);
   if(b>0) sum += ramp_sumaccum_gauss(SUM_GAUSS_ARGS(vector,max(0,a),b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}

vector ramp_slope_hold_gauss( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_GAUSS_HEAD(vector);
   if(a<0)        sum += ramp_sumslope_gauss(0,lv,lm,-0.5,(min(b,0)-x)*iw,x,w);
   if(a<1 && b>0) 
      sum += ramp_sum_gauss(nl,VZERO,max(a,0),min(b,1),x,w,iw,bases,keys,vals);
   if(b>1)        sum += ramp_sumconst_gauss(rv,(max(a,1)-x)*iw,0.5);
   return sum*RAMP_GAUSS_INORM;
}
vector ramp_slope_slope_gauss( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_GAUSS_HEAD(vector);
   if(a<0)        sum += ramp_sumslope_gauss(0,lv,lm,-0.5,(min(b,0)-x)*iw,x,w);
   if(a<1 && b>0) 
      sum += ramp_sum_gauss(nl,VZERO,max(a,0),min(b,1),x,w,iw,bases,keys,vals);
   if(b>1)        sum += ramp_sumslope_gauss(1,rv,rm,(max(a,1)-x)*iw,0.5,x,w);
   return sum*RAMP_GAUSS_INORM;
}
vector ramp_slope_cycle_gauss( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_GAUSS_HEAD(vector);
   if(a<0)        sum += ramp_sumslope_gauss(0,lv,lm,-0.5,(min(b,0)-x)*iw,x,w);
   if(b>0) sum += ramp_sumcycle_gauss(SUM_GAUSS_ARGS(vector,max(0,a),b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}
vector ramp_slope_mirror_gauss( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_GAUSS_HEAD(vector);
   if(a<0)        sum += ramp_sumslope_gauss(0,lv,lm,-0.5,(min(b,0)-x)*iw,x,w);
   if(b>0) sum += ramp_summirror_gauss(SUM_GAUSS_ARGS(vector,max(0,a),b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}
vector ramp_slope_accum_gauss( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_GAUSS_HEAD(vector);
   if(a<0)        sum += ramp_sumslope_gauss(0,lv,lm,-0.5,(min(b,0)-x)*iw,x,w);
   if(b>0) sum += ramp_sumaccum_gauss(SUM_GAUSS_ARGS(vector,max(0,a),b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}

vector ramp_cycle_hold_gauss( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_GAUSS_HEAD(vector);
   if(a<1) sum += ramp_sumcycle_gauss(SUM_GAUSS_ARGS(vector,a,min(1,b),x,w) );
   if(b>1)        sum += ramp_sumconst_gauss(rv,(max(a,1)-x)*iw,0.5);
   return sum*RAMP_GAUSS_INORM;
}
vector ramp_cycle_slope_gauss( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_GAUSS_HEAD(vector);
   if(a<1) sum += ramp_sumcycle_gauss(SUM_GAUSS_ARGS(vector,a,min(1,b),x,w) );
   if(b>1)        sum += ramp_sumslope_gauss(1,rv,rm,(max(a,1)-x)*iw,0.5,x,w);
   return sum*RAMP_GAUSS_INORM;
}
vector ramp_cycle_cycle_gauss( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_GAUSS_HEAD(vector);
   sum = ramp_sumcycle_gauss(SUM_GAUSS_ARGS(vector,a,b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}
vector ramp_cycle_mirror_gauss( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_GAUSS_HEAD(vector);
   if(a<1) sum += ramp_sumcycle_gauss(SUM_GAUSS_ARGS(vector,a,min(1,b),x,w) );
   if(b>1) sum += ramp_summirror_gauss(SUM_GAUSS_ARGS(vector,max(1,a),b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}
vector ramp_cycle_accum_gauss( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_GAUSS_HEAD(vector);
   if(a<1) sum += ramp_sumcycle_gauss(SUM_GAUSS_ARGS(vector,a,min(1,b),x,w) );
   if(b>1) sum += ramp_sumaccum_gauss(SUM_GAUSS_ARGS(vector,max(1,a),b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}

vector ramp_mirror_hold_gauss( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_GAUSS_HEAD(vector);
   if(a<1) sum += ramp_summirror_gauss(SUM_GAUSS_ARGS(vector,a,min(1,b),x,w) );
   if(b>1) sum += ramp_sumconst_gauss(rv,(max(a,1)-x)*iw,0.5);
   return sum*RAMP_GAUSS_INORM;
}
vector ramp_mirror_slope_gauss( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_GAUSS_HEAD(vector);
   if(a<1) sum += ramp_summirror_gauss(SUM_GAUSS_ARGS(vector,a,min(1,b),x,w) );
   if(b>1) sum += ramp_sumslope_gauss(1,rv,rm,(max(a,1)-x)*iw,0.5,x,w);
   return sum*RAMP_GAUSS_INORM;
}
vector ramp_mirror_cycle_gauss( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_GAUSS_HEAD(vector);
   if(a<1) sum += ramp_summirror_gauss(SUM_GAUSS_ARGS(vector,a,min(1,b),x,w) );
   if(b>1) sum += ramp_sumcycle_gauss(SUM_GAUSS_ARGS(vector,max(1,a),b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}
vector ramp_mirror_mirror_gauss( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_GAUSS_HEAD(vector);
   sum += ramp_summirror_gauss(SUM_GAUSS_ARGS(vector,a,b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}
vector ramp_mirror_accum_gauss( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_GAUSS_HEAD(vector);
   if(a<1) sum += ramp_summirror_gauss(SUM_GAUSS_ARGS(vector,a,min(1,b),x,w) );
   if(b>1) sum += ramp_sumaccum_gauss(SUM_GAUSS_ARGS(vector,max(1,a),b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}

vector ramp_accum_hold_gauss( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_GAUSS_HEAD(vector);
   if(a<1) sum += ramp_sumaccum_gauss(SUM_GAUSS_ARGS(vector,a,min(1,b),x,w) );
   if(b>1) sum += ramp_sumconst_gauss(rv,(max(a,1)-x)*iw,0.5);
   return sum*RAMP_GAUSS_INORM;
}
vector ramp_accum_slope_gauss( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_GAUSS_HEAD(vector);
   if(a<1) sum += ramp_sumaccum_gauss(SUM_GAUSS_ARGS(vector,a,min(1,b),x,w) );
   if(b>1) sum += ramp_sumslope_gauss(1,rv,rm,(max(a,1)-x)*iw,0.5,x,w);
   return sum*RAMP_GAUSS_INORM;
}
vector ramp_accum_cycle_gauss( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_GAUSS_HEAD(vector);
   if(a<1) sum += ramp_sumaccum_gauss(SUM_GAUSS_ARGS(vector,a,min(1,b),x,w) );
   if(b>1) sum += ramp_sumcycle_gauss(SUM_GAUSS_ARGS(vector,max(1,a),b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}
vector ramp_accum_mirror_gauss( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_GAUSS_HEAD(vector);
   if(a<1) sum += ramp_sumaccum_gauss(SUM_GAUSS_ARGS(vector,a,min(1,b),x,w) );
   if(b>1) sum += ramp_summirror_gauss(SUM_GAUSS_ARGS(vector,max(1,a),b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}
vector ramp_accum_accum_gauss( RAMP_FILTERED_PARMS(vector) ) {
   RAMP_GAUSS_HEAD(vector);
   sum += ramp_sumaccum_gauss(SUM_GAUSS_ARGS(vector,a,b,x,w) );
   return sum*RAMP_GAUSS_INORM;
}


// Universal functions which can be used in VOPs without requiring any
// recompilation due to parm changes.
#define EXTRBLOCK(EXTL, EXTR, FILTER) \
(extr == 'EXTR') \
{ \
    return ramp_##EXTL##_##EXTR##_##FILTER(RAMP_FILTERED_ARGS); \
}

#define EXTLBLOCK(EXTL, FILTER) \
if(extl == 'EXTL') \
{ \
    if	    EXTRBLOCK(EXTL,   hold,   FILTER) \
    else if EXTRBLOCK(EXTL,   slope,  FILTER) \
    else if EXTRBLOCK(EXTL,   cycle,  FILTER) \
    else if EXTRBLOCK(EXTL,   mirror, FILTER) \
    else if EXTRBLOCK(EXTL,   accum,  FILTER) \
}

#define RAMPRTYPEBLOCK(RTYPE) \
RTYPE ramp_lookup(string extl, extr, filter; RAMP_FILTERED_PARMS(RTYPE)) \
{ \
    if(filter == "point") \
    { \
	EXTLBLOCK(hold, point) \
	EXTLBLOCK(slope, point) \
	EXTLBLOCK(cycle, point) \
	EXTLBLOCK(mirror, point) \
	EXTLBLOCK(accum, point) \
    } \
    else if(filter == "box") \
    { \
	EXTLBLOCK(hold, box) \
	EXTLBLOCK(slope, box) \
	EXTLBLOCK(cycle, box) \
	EXTLBLOCK(mirror, box) \
	EXTLBLOCK(accum, box) \
    } \
    else if(filter == "gauss") \
    { \
	EXTLBLOCK(hold, gauss) \
	EXTLBLOCK(slope, gauss) \
	EXTLBLOCK(cycle, gauss) \
	EXTLBLOCK(mirror, gauss) \
	EXTLBLOCK(accum, gauss) \
    } \
    return 0.0; \
}

RAMPRTYPEBLOCK(float)
RAMPRTYPEBLOCK(vector)

#endif // End pyro_aaramp_h
