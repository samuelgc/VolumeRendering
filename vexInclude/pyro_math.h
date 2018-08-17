/****************************************************************************
   Description:  Basic math support for the pyro tools

*****************************************************************************/

#ifndef pyro_math_h__GUARD
#define pyro_math_h__GUARD
   

#include <pyro_types.h>


// -------------------------------------------------------------------------------
// Assorted constants, renamed from M_ to avoid clashes
// -------------------------------------------------------------------------------
#define C_PI               3.14159265358979323846  // pi 
#define C_PI2              6.28318530717958647692  // pi*2
#define C_PI_2             1.57079632679489661923  // pi/2
#define C_PI4              12.5663706143591729539  // pi*4
#define C_PI_4             0.78539816339744830962  // pi/4
#define C_SQRT2            1.41421356237309504880  // sqrt(2) 
#define C_1_SQRT2          0.70710678118654752440  // 1/sqrt(2)
#define C_SQRT3            1.73205080756887729353  // sqrt(3)
#define C_1_SQRT3          0.57735026918962576451  // 1/sqrt(3)
#define C_SQRTPI           1.7724538509055160273   // sqrt(PI)
#define C_SQRT2PI          2.50662827463100050242  // sqrt(2*PI)
#define C_SQRTPI_2         1.25331413731550025121  // sqrt(PI/2)
#define C_1_3              0.3333333333333333333   // 1/3
#define C_LOG2             0.69314718055994530942  // log(2)
#define C_E                2.71828182845904523536  // e
#define C_E2               5.43656365691809047072  // 2*e

#define C_FLOATEPSILON     1e-5

#define VONE               {1,1,1}
#define VTWO               {2,2,2}
#define VZERO              {0,0,0}

// -------------------------------------------------------------------------------
// Some low level support extending vex/include/math.h
// -------------------------------------------------------------------------------
void sort2(export float a,b) {
   if(a>b) {
      float tmp = a;
      a = b;
      b = tmp;
   }
}

float modulo(float a,b) { 
   float r=a%b; return r<0 ? r+b : r; 
}
vector modulo(vector a,b) { 
   return set(modulo(a.x,b.x),
              modulo(a.y,b.y),
              modulo(a.z,b.z)); 
}
vector modulo(vector a; float b) { 
   return set(modulo(a.x,b),
              modulo(a.y,b),
              modulo(a.z,b)); 
}
vector4 modulo(vector4 a,b) { 
   return set(modulo(a.x,b.x),
              modulo(a.y,b.y),
              modulo(a.z,b.z), 
              modulo(a.w,b.w)); 
}
vector4 modulo(vector4 a; float b) { 
   return set(modulo(a.x,b),
              modulo(a.y,b),
              modulo(a.z,b), 
              modulo(a.w,b)); 
}


int equalzero(float val) {
   return (val >= -C_FLOATEPSILON && val<= C_FLOATEPSILON);
}
int equalzero(float val,eps) {
   return (val >= -eps && val<= eps);
}
int equalzero(vector val) {
   return (equalzero(val.x) && equalzero(val.y) && equalzero(val.z));
}

int exceptzero(int x) {
   return x == 0 ? 1: x;
}
float exceptzero(float x) {
   return equalzero(x,C_FLOATEPSILON) ? C_FLOATEPSILON*sign(x) : x;
}
float exceptzero(float x,eps) {
   return equalzero(x,eps) ? eps*sign(x) : x;
}
vector exceptzero(vector x) {
   return set(exceptzero(x.x), exceptzero(x.y), exceptzero(x.z));
}
vector exceptzero(vector x; float eps) {
   return set(exceptzero(x.x,eps), exceptzero(x.y,eps), exceptzero(x.z,eps));
}
vector exceptzero(vector x, eps) {
   return set(exceptzero(x.x,eps.x), exceptzero(x.y,eps.y), exceptzero(x.z,eps.z));
}
vector4 exceptzero(vector4 x) {
   return set(exceptzero(x.x), exceptzero(x.y), exceptzero(x.z), exceptzero(x.w));
}
vector4 exceptzero(vector4 x; float eps) {
   return set(exceptzero(x.x,eps), exceptzero(x.y,eps), 
              exceptzero(x.z,eps), exceptzero(x.w,eps));
}
vector4 exceptzero(vector4 x, eps) {
   return set(exceptzero(x.x,eps.x), exceptzero(x.y,eps.y), 
              exceptzero(x.z,eps.z), exceptzero(x.w,eps.w));
}

float fabs_integral(float t) {
   return sign(t) * 0.5 * t*t;
}
float fabs(float x, dx) {
   float x0=x-0.5*dx;
   float x1=x0+dx;
   return (fabs_integral(x1)-fabs_integral(x0))/dx;
}
#define fabsf fabs
vector fabsv(vector x; float dx) {
   return set( fabs(x[0],dx), fabs(x[1],dx), fabs(x[2],dx) );
}
vector4 fabsv4(vector4 x; float dx) {
   return set( fabs(x[0],dx), fabs(x[1],dx), fabs(x[2],dx), fabs(x[3],dx) );
}



vector lerpv(vector a,b,bias) {
   return set(lerp(a.x,b.x,bias.x),
              lerp(a.y,b.y,bias.y),
              lerp(a.z,b.z,bias.z));
}

vector smoothv(vector a,b,bias) {
   return set(smooth(a.x,b.x,bias.x),
              smooth(a.y,b.y,bias.y),
              smooth(a.z,b.z,bias.z));
}


// Extended max()
vector max(vector a; float b) {
   return set ( max(a.x,b), max(a.y,b), max(a.z,b) );
}
vector max(float a; vector b) {
   return max(b,a);
}
vector4 maxv(vector4 a; float b) {
   return set( max(a.x,b), max(a.y,b), max(a.z,b), max(a.w,b) );
}
vector4 maxf(float a; vector4 b) {
   return maxv(b,a);
}
int[] max(int a[]; int b[]) {
   int out[] = array();
   resize(out,min(len(a),len(b)));
   int i, l = len(out);
   for(i=0;i<l;i++) out[i] = max(a[i],b[i]);
   return out;
}
int[] max(int a[]; int b) {
   int out[] = a;
   int i, l = len(out);
   for(i=0;i<l;i++) out[i] = max(a[i],b);
   return out;
}
int[] maxiia(int a; int b[]) {
   return max(b,a);
}
float[] max(float a[]; float b[]) {
   float out[] = array(); 
   resize(out,min(len(a),len(b)));
   int i, l = len(out);
   for(i=0;i<l;i++) out[i] = max(a[i],b[i]);
   return out;
}
float[] max(float a[]; float b) {
   float out[] = a;
   int i, l = len(out);
   for(i=0;i<l;i++) out[i] = max(a[i],b);
   return out;
}
float[] max(float a; float b[]) {
   return max(b,a);
}
vector[] max(vector a[]; vector b[]) {
   vector out[] = array(); 
   resize(out,min(len(a),len(b)));
   int i, l = len(out);
   for(i=0;i<l;i++) out[i] = max(a[i],b[i]);
   return out;
}
vector[] max(vector a[]; vector b) {
   vector out[] = a;
   int i, l = len(out);
   for(i=0;i<l;i++) out[i] = max(a[i],b);
   return out;
}
vector[] max(vector a; vector b[]) {
   return max(b,a);
}
vector[] max(vector a[]; float b) {
   vector out[] = a;
   int i, l = len(out);
   for(i=0;i<l;i++) out[i] = max(a[i],b);
   return out;
}
vector[] max(float a; vector b[]) {
   return max(b,a);
}
vector4[] max(vector4 a[]; vector4 b[]) {
   vector4 out[] = array(); 
   resize(out,min(len(a),len(b)));
   int i, l = len(out);
   for(i=0;i<l;i++) out[i] = max(a[i],b[i]);
   return out;
}
vector4[] max(vector4 a[]; vector4 b) {
   vector4 out[] = a;
   int i, l = len(out);
   for(i=0;i<l;i++) out[i] = max(a[i],b);
   return out;
}
vector4[] max(vector4 a; vector4 b[]) {
   return max(b,a);
}
vector4[] max(vector4 a[]; float b) {
   vector4 out[] = a;
   int i, l = len(out);
   for(i=0;i<l;i++) out[i] = max(a[i],b);
   return out;
}
vector4[] max(float a; vector4 b[]) {
   return max(b,a);
}




// Extended min()
vector min(vector a; float b) {
   return set ( min(a.x,b), min(a.y,b), min(a.z,b) );
}
vector min(float a; vector b) {
   return min(b,a);
}
vector4 min(vector4 a; float b) {
   return set ( min(a.x,b), min(a.y,b), min(a.z,b), min(a.w,b) );
}
vector4 min(float a; vector4 b) {
   return min(b,a);
}
int[] min(int a[]; int b[]) {
   int out[] = array();
   resize(out,min(len(a),len(b)));
   int i, l = len(out);
   for(i=0;i<l;i++) out[i] = min(a[i],b[i]);
   return out;
}
int[] min(int a[]; int b) {
   int out[] = a;
   int i, l = len(out);
   for(i=0;i<l;i++) out[i] = min(a[i],b);
   return out;
}
int[] min(int a; int b[]) {
   return min(b,a);
}
float[] min(float a[]; float b[]) {
   float out[] = array(); 
   resize(out,min(len(a),len(b)));
   int i, l = len(out);
   for(i=0;i<l;i++) out[i] = min(a[i],b[i]);
   return out;
}
float[] min(float a[]; float b) {
   float out[] = a;
   int i, l = len(out);
   for(i=0;i<l;i++) out[i] = min(a[i],b);
   return out;
}
float[] min(float a; float b[]) {
   return min(b,a);
}
vector[] min(vector a[]; vector b[]) {
   vector out[] = array(); 
   resize(out,min(len(a),len(b)));
   int i, l = len(out);
   for(i=0;i<l;i++) out[i] = min(a[i],b[i]);
   return out;
}
vector[] min(vector a[]; vector b) {
   vector out[] = a;
   int i, l = len(out);
   for(i=0;i<l;i++) out[i] = min(a[i],b);
   return out;
}
vector[] min(vector a; vector b[]) {
   return min(b,a);
}
vector[] min(vector a[]; float b) {
   vector out[] = a;
   int i, l = len(out);
   for(i=0;i<l;i++) out[i] = min(a[i],b);
   return out;
}
vector[] min(float a; vector b[]) {
   return min(b,a);
}
vector4[] min(vector4 a[]; vector4 b[]) {
   vector4 out[] = array(); 
   resize(out,min(len(a),len(b)));
   int i, l = len(out);
   for(i=0;i<l;i++) out[i] = min(a[i],b[i]);
   return out;
}
vector4[] min(vector4 a[]; vector4 b) {
   vector4 out[] = a;
   int i, l = len(out);
   for(i=0;i<l;i++) out[i] = min(a[i],b);
   return out;
}
vector4[] min(vector4 a; vector4 b[]) {
   return min(b,a);
}
vector4[] min(vector4 a[]; float b) {
   vector4 out[] = a;
   int i, l = len(out);
   for(i=0;i<l;i++) out[i] = min(a[i],b);
   return out;
}
vector4[] min(float a; vector4 b[]) {
   return min(b,a);
}



// extended pow
vector pow(float b; vector e) {
   return set(pow(b,e.x), pow(b,e.y), pow(b,e.z));
}
vector pow(vector b,e) {
   return set(pow(b.x,e.x), pow(b.y,e.y), pow(b.z,e.z));
}
vector4 pow(float b; vector4 e) {
   return set(pow(b,e.x), pow(b,e.y), pow(b,e.z), pow(b,e.w));
}
vector4 pow(vector4 b,e) {
   return set(pow(b.x,e.x), pow(b.y,e.y), pow(b.z,e.z), pow(b.w,e.w));
}
#define ARRAYAPPLYPOW(outtype,btype,etype) \
outtype[] pow(btype b[]; etype e)   {  \
   int len = len(b); \
   outtype out[]; \
   if(len<=0) return out; \
   resize(out,len); \
   int i; for(i=0;i<len;i++) out[i] = pow(b[i],e); \
   return out; \
}
ARRAYAPPLYPOW(float,float,float)
ARRAYAPPLYPOW(vector,vector,float)
ARRAYAPPLYPOW(vector4,vector4,float)
ARRAYAPPLYPOW(vector,vector,vector)
ARRAYAPPLYPOW(vector4,vector4,vector4)
ARRAYAPPLYPOW(vector,float,vector)
ARRAYAPPLYPOW(vector4,float,vector4)
#undef ARRAYAPPLYPOW

#define ARRAYAPPLYPROC(type,fin,initval,fout) \
type fin(type x[])   {  \
   int len = len(x); \
   if(len<=0) return 0; \
   type out = type(initval); \
   int i; for(i=0;i<len;i++) out = (type)fin(out,x[i]); \
   float fi = (float)i; \
   out = (type)(fout); \
   return out; \
}
#define ARRAYAPPLY(fin,initval,fout) \
   ARRAYAPPLYPROC(int     , fin , initval, fout) \
   ARRAYAPPLYPROC(vector  , fin , initval, fout) \
   ARRAYAPPLYPROC(vector4 , fin , initval, fout)


int     sum(int a,b)     { return a+b; }
float   sum(float a,b)   { return a+b; }
vector  sum(vector a,b)  { return a+b; }
vector4 sum(vector4 a,b) { return a+b; }

#undef ARRAYAPPLY
#undef ARRAYAPPLYPROC


// component product
#define ARRAYPROD(type) \
type prod(type x[])   {  \
   int len = len(x); \
   if(len<=0) return (type)0; \
   type out = (type)1; \
   int i; for(i=0;i<len;i++) out *= x[i]; \
   return out; \
}
int   prod(int x)     { return x; }
float prod(float x)   { return x; }
float prod(vector x)  { return x.x*x.y*x.z; }
float prod(vector4 x) { return x.x*x.y*x.z*x.w; }
ARRAYPROD(int)
ARRAYPROD(vector)
ARRAYPROD(vector4)
#undef ARRAYPROD

// extended log (log to a given base)
float  log(float x, base) {
   return log(max(1e-37,x)) / log(max(1e-37,base));
}
vector log(vector x; float base) {
   return log(max(1e-37,x)) / log(max(1e-37,base));
}
vector log(vector x, base) {
   return log(max(1e-37,x)) / log(max(1e-37,base));
}
vector4 log(vector4 x; float base) {
   return log(max(1e-37,x)) / log(max(1e-37,base));
}
vector4 log(vector4 x, base) {
   return log(max(1e-37,x)) / log(max(1e-37,base));
}

// Ken Perlin's Bias function
// The bias parameter 'b' must be in (0,1)
float bias(float x, b) {
   if(b==0.5) return x;
   float bb = clamp(b,1e-6,1-1e-6);
   return pow(x,-log(bb)/C_LOG2);
}
vector bias(vector x; float b) {
   if(b==0.5) return x;
   float bb = clamp(b,1e-6,1-1e-6);
   return pow(x,-log(bb)/C_LOG2);
}
vector bias(vector x, b) {
   if(b==(vector)0.5) return x;
   return set(bias(x.x,b.x),bias(x.y,b.y),bias(x.z,b.z));
}




// Ken Perlin's Gain function -- really only valid for x in [0,1]
// The gain parameter 'g' must be in (0,1)
float gain(float x,g) {
   if(g==0.5) return x;
   float gc = 1-g;
   return (x<.5) ? bias(2*x,gc)*0.5 : 1 - bias(2-2*x,gc)*0.5;
}
vector gain(vector x; float g) {
   if(g==0.5) return x;
   float  gc = 1-g;
   vector x2 = x*2;
   return set(
         (x.x<.5) ? bias(x2.x,gc)*0.5 : 1 - bias(2-x2.x,gc)*0.5,
         (x.y<.5) ? bias(x2.y,gc)*0.5 : 1 - bias(2-x2.y,gc)*0.5,
         (x.z<.5) ? bias(x2.z,gc)*0.5 : 1 - bias(2-x2.z,gc)*0.5 );
}
vector gain(vector x,g) {
   if(g==(vector)0.5) return x;
   vector gc = VONE-g;
   vector x2 = x*2;
   return set(
         (x.x<.5) ? bias(x2.x,gc.x)*0.5 : 1 - bias(2-x2.x,gc.x)*0.5,
         (x.y<.5) ? bias(x2.y,gc.y)*0.5 : 1 - bias(2-x2.y,gc.y)*0.5,
         (x.z<.5) ? bias(x2.z,gc.z)*0.5 : 1 - bias(2-x2.z,gc.z)*0.5 );
}



float ndxmirror(float x) {
   return 1.0 - 2.0*abs(0.5-modulo(x*.5,1));
}
vector ndxmirror(vector x) {
   return set(ndxmirror(x.x),ndxmirror(x.y),ndxmirror(x.z));
}
vector4 ndxmirror(vector4 x) {
   return set(ndxmirror(x.x),ndxmirror(x.y),ndxmirror(x.z),ndxmirror(x.w));
}



// erf() and erfc() have now been added to VEX, so the following is not
// needed any more. I'm keeping it around a little longer just in case.
//------------------------------------------------------------------------------

#if 0
// A few different approaches to approximating the erf() function (used
// elsewhere in pyro to integrate or convolve with Gaussians).
float erf_hmf(float x) {
   if(x==0) return 0;

   float sx = sign(x);
   float ax = x*sx;
   
   if(ax<1e-10) 
      return sx*(ax*1.125 + ax*0.003379167095512573896158903121545171688);

   float a1 =  0.254829592;
   float a2 = -0.284496736;
   float a3 =  1.421413741;
   float a4 = -1.453152027;
   float a5 =  1.061405429;
   float p  =  0.3275911;
   
   float k = 1.0/(1.0 + p*ax);
   float y = 1.0 - ((((a5*k + a4)*k + a3)*k + a2)*k + a1)*k*exp(-ax*ax);
   
   return sx*y;
}

float erf_wp(float x) {
   if(x==0) return 0;

   float sx = sign(x);
   float ax = x*sx;
   
   if(ax<1e-10) 
      return sx*(ax*1.125 + ax*0.003379167095512573896158903121545171688);

   float ax2 = ax*ax;
   float a   = 0.140012288686666606004249;
   return sx*sqrt(1.0-exp(-ax2*(4.0+a*C_PI*ax2)/(C_PI+a*C_PI*ax2)));
}

float erf_a(float x) {
   float sx = sign(x);
   float ax = x*sx;

   if(ax<1e-10) 
      return sx*(ax*1.125 + ax*0.003379167095512573896158903121545171688);

   float k  = 1/(1+0.5*ax);
   return sx * (
      1-k*exp(-ax*ax - 1.26551223 +
                 k * ( 1.00002368 +
                 k * ( 0.37409196 + 
                 k * ( 0.09678418 + 
                 k * (-0.18628806 + 
                 k * ( 0.27886807 + 
                 k * (-1.13520398 + 
                 k * ( 1.48851587 + 
                 k * (-0.82215223 + 
                 k * ( 0.17087277))))))))))
   );
}



// As a lookup table. erf(>4.5) is forced to 1.0.
float erf_lookup(float x) {
   float e[] = {
      0.000000000000, 0.016981013038, 0.033954335454, 0.050912287073, 0.067847208591, 
      0.084751471954, 0.101617490661, 0.118437729977, 0.135204717029, 0.151911050763, 
      0.168549411736, 0.185112571731, 0.201593403161, 0.217984888256, 0.234280128000, 
      0.250472350809, 0.266554920923, 0.282521346506, 0.298365287429, 0.314080562718, 
      0.329661157666, 0.345101230581, 0.360395119167, 0.375537346524, 0.390522626755, 
      0.405345870175, 0.420002188112, 0.434486897296, 0.448795523821, 0.462923806694, 
      0.476867700952, 0.490623380349, 0.504187239619, 0.517555896308, 0.530726192182, 
      0.543695194205, 0.556460195110, 0.569018713542, 0.581368493802, 0.593507505182, 
      0.605433940919, 0.617146216752, 0.628642969107, 0.639923052928, 0.650985539147, 
      0.661829711811, 0.672455064894, 0.682861298785, 0.693048316482, 0.703016219498, 
      0.712765303501, 0.722296053696, 0.731609139974, 0.740705411830, 0.749585893088, 
      0.758251776426, 0.766704417731, 0.774945330304, 0.782976178919, 0.790798773763, 
      0.798415064269, 0.805827132865, 0.813037188641, 0.820047560967, 0.826860693070, 
      0.833479135579, 0.839905540065, 0.846142652585, 0.852193307245, 0.858060419792, 
      0.863746981255, 0.869256051645, 0.874590753722, 0.879754266849, 0.884749820940, 
      0.889580690505, 0.894250188817, 0.898761662201, 0.903118484449, 0.907324051383, 
      0.911381775559, 0.915295081131, 0.919067398865, 0.922702161334, 0.926202798264, 
      0.929572732072, 0.932815373574, 0.935934117867, 0.938932340411, 0.941813393279, 
      0.944580601606, 0.947237260214, 0.949786630434, 0.952231937105, 0.954576365769, 
      0.956823060042, 0.958975119167, 0.961035595761, 0.963007493718, 0.964893766311, 
      0.966697314445, 0.968420985091, 0.970067569881, 0.971639803863, 0.973140364410, 
      0.974571870281, 0.975936880836, 0.977237895377, 0.978477352642, 0.979657630418, 
      0.980781045286, 0.981849852485, 0.982866245891, 0.983832358106, 0.984750260648, 
      0.985621964244, 0.986449419213, 0.987234515935, 0.987979085407, 0.988684899870, 
      0.989353673507, 0.989987063217, 0.990586669433, 0.991154037012, 0.991690656166, 
      0.992197963440, 0.992677342738, 0.993130126378, 0.993557596182, 0.993960984600, 
      0.994341475848, 0.994700207076, 0.995038269548, 0.995356709834, 0.995656531019, 
      0.995938693906, 0.996204118236, 0.996453683902, 0.996688232157, 0.996908566830, 
      0.997115455518, 0.997309630786, 0.997491791344, 0.997662603213, 0.997822700882, 
      0.997972688442, 0.998113140702, 0.998244604294, 0.998367598746, 0.998482617540, 
      0.998590129151, 0.998690578053, 0.998784385707, 0.998871951530, 0.998953653824, 
      0.999029850695, 0.999100880935, 0.999167064888, 0.999228705279, 0.999286088025, 
      0.999339483018, 0.999389144882, 0.999435313701, 0.999478215723, 0.999518064046, 
      0.999555059263, 0.999589390099, 0.999621234014, 0.999650757783, 0.999678118054, 
      0.999703461884, 0.999726927253, 0.999748643551, 0.999768732048, 0.999787306345, 
      0.999804472798, 0.999820330931, 0.999834973819, 0.999848488467, 0.999860956158, 
      0.999872452787, 0.999883049186, 0.999892811424, 0.999901801092, 0.999910075583, 
      0.999917688342, 0.999924689119, 0.999931124194, 0.999937036600, 0.999942466329, 
      0.999947450528, 0.999952023681, 0.999956217788, 0.999960062524, 0.999963585394, 
      0.999966811883, 0.999969765588, 0.999972468345, 0.999974940355, 0.999977200295, 
      0.999979265422, 0.999981151674, 0.999982873765, 0.999984445270, 0.999985878708, 
      0.999987185615, 0.999988376622, 0.999989461515, 0.999990449301, 0.999991348265, 
      0.999992166022, 0.999992909573, 0.999993585343, 0.999994199235, 0.999994756661, 
      0.999995262586, 0.999995721560, 0.999996137752, 0.999996514978, 0.999996856732, 
      0.999997166210, 0.999997446332, 0.999997699768, 0.999997928958, 0.999998136125, 
      0.999998323303, 0.999998492343, 0.999998644933, 0.999998782612, 0.999998906780, 
      0.999999018714, 0.999999119572, 0.999999210410, 0.999999292185, 0.999999365769, 
      0.999999431953, 0.999999491452, 0.999999544919, 0.999999592944, 0.999999636059, 
      0.999999674751, 0.999999709457, 0.999999740573, 0.999999768459, 0.999999793437, 
      0.999999815802, 0.999999835818, 0.999999853723, 0.999999869733, 0.999999884041, 
      0.999999896823, 0.999999908236, 0.999999918423, 0.999999927511, 0.999999935615, 
      0.999999942838, 0.999999949273, 0.999999955003, 0.999999960104, 0.999999964642, 
      0.999999968678, 0.999999972265, 0.999999975452, 0.999999978283, 0.999999980795, 
      0.999999983025, 0.999999985002, 0.999999986755, 0.999999988308, 0.999999989683, 
      0.999999990901, 0.999999991978, 0.999999992931, 0.999999993774, 0.999999994519, 
      0.999999995176, 0.999999995757, 0.999999996269, 0.999999996721, 0.999999997120, 
      0.999999997471, 0.999999997780, 0.999999998053, 0.999999998292, 0.999999998503, 
      0.999999998689, 0.999999998851, 0.999999998995, 0.999999999120, 0.999999999231, 
      0.999999999327, 0.999999999412, 0.999999999487, 0.999999999552, 0.999999999609, 
      0.999999999659, 0.999999999703, 0.999999999741, 0.999999999774, 1.000000000000
   };
   float sx  = sign(x);
   float ax  = x*sx; // abs(x)
   if(ax>=4.5) return sx;

   //return spline("linear",ax/4.5,e)*sx;

   // note: this is ~2X faster than spline("linear",ax/4.5,e)*sx
   float xf = ax*66.444444444444; // == ax*(len(e)-1)/4.5
   int   xi = floor(xf);
   xf -= xi;
   return sx*lerp(e[xi],e[xi+1],xf);
}


// Pick one of the erf() implementations
#define erf erf_hmf


float erfc(float x) {
   return 1.0 - erf(x);
}
#endif


#endif // End pyro_math_h
