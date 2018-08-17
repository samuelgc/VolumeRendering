/****************************************************************************
   Description:  The shading module for pyro2

*****************************************************************************/

#ifndef pyro_shade_h__GUARD
#define pyro_shade_h__GUARD

#include <pyro_utils.h>
#include <cie_cmf.h>



//--------------------------------------------------------------------------------
// XYZ and RGB
//--------------------------------------------------------------------------------

#define WP_D65        {0.3127, 0.3290, 0.3583} // White Point: D65    z = 0.3583
#define LUMA_Rec709   {0.2126, 0.7152, 0.0722} // Luma: HDTV          z = 0.0722
#define LUMA_Rec601   {0.299 , 0.587 , 0.114}  // Luma: SDTV          z = 0.1140
#define LUMA_Haeberli {0.3086, 0.6094, 0.0820} // Luma: Paul Haeberli z = 0.0820

#define LUMA_Default   LUMA_Haeberli           // The default luma metric
                                               // ..just coz it matches COPs..

// Gamut {x,y} for sRGB (ITU-R BT.709 Primaries, WP=D65)
#define GAMUT_R   {0.6400, 0.3300, 0}        // z = (1-x-y) = 0.0300
#define GAMUT_G   {0.3000, 0.6000, 0}        // z = (1-x-y) = 0.1000
#define GAMUT_B   {0.1500, 0.0600, 0}        // z = (1-x-y) = 0.7900

// Barycentric and line-intersect constants for the sRGB gamut
#define GAMUT_RB  {-0.49, -0.27, 0}          // = (GAMUT_B-GAMUT_R)
#define GAMUT_RG  {-0.34,  0.27, 0}          // = (GAMUT_G-GAMUT_R)
#define GAMUT_GB  {-0.15, -0.54, 0}          // = (GAMUT_B-GAMUT_G)
#define GAMUT_LRB 0.559464                   // = length(GAMUT_RB)
#define GAMUT_LRG 0.434166                   // = length(GAMUT_RG)
#define GAMUT_LGB 0.560446                   // = length(GAMUT_GB)
#define GAMUT_NRB {-0.875838, -0.482605, 0}  // = normalize(GAMUT_RB)
#define GAMUT_NRG {-0.783111,  0.621882, 0}  // = normalize(GAMUT_RG)
#define GAMUT_NGB {-0.267644, -0.963518, 0}  // = normalize(GAMUT_GB)
#define GAMUT_RW  {-0.3273, -0.001, 0}       // = (WP_D65-GAMUT_R)
#define GAMUT_GW  { 0.0127, -0.271, 0}       // = (WP_D65-GAMUT_G)
#define GAMUT_BW  { 0.1627,  0.269, 0}       // = (WP_D65-GAMUT_B)
#define GAMUT_B00 0.313                      // = V0 dot self
#define GAMUT_B01 0.0937                     // = V0 dot V1
#define GAMUT_B11 0.1885                     // = V1 dot self
#define GAMUT_BK  19.9121                    // = 1 / (B00*B11 - B01*B01)

// video luminance (excluding gamma)
float luma(vector linrgb) {
   return dot(linrgb,LUMA_Default);
}

// perceived lightness
float lightness(vector linrgb) {
   return (ctransform("cspace:Lab",linrgb))[0];
}

// XYZ -> xyY
vector xyztoxyy(vector xyz) {
   float d = xyz.x + xyz.y + xyz.z;
   if(d==0) return 0;
   return set(xyz.x/d, xyz.y/d, xyz.y);
}

// XYZ -> xyY (exports norm: L1 metric)
vector xyztoxyy(vector xyz; export float d) {
   d = xyz.x + xyz.y + xyz.z;
   if(d==0) return 0;
   return set(xyz.x/d, xyz.y/d, xyz.y);
}

// xyY -> XYZ
vector xyytoxyz(vector xyY) {
   //if(xyY.y<=0) return 0;
   float k = xyY.z / xyY.y;
   return set( k * xyY.x, xyY.z, k * (1.0 - xyY.x - xyY.y) );
}

// XYZ -> Chromaticity[3] (xyz)
vector xyztochroma(vector xyz) {
   float d = xyz.x + xyz.y + xyz.z;
   if(d==0) return 0;
   return xyz / d;
}

// xyY -> Linear RGB
vector xyytorgb(vector xyy) {
   return xyztorgb(xyytoxyz(xyy));
}

// Linear RGB -> xyY
vector rgbtoxyy(vector rgb) {
   return xyztoxyy(rgbtoxyz(rgb));
}

// Linear RGB -> xyY (exports norm: L1 metric)
vector rgbtoxyy(vector rgb; export float d) {
   return xyztoxyy(rgbtoxyz(rgb),d);
}



//--------------------------------------------------------------------------------
// Gamut Mapping
//--------------------------------------------------------------------------------

// Barycentric coords of chromaticity {x,y} relative to sRGB gamut
void barysrgb(vector xy; export float bu; export float bv) {
   vector v2  = xy - GAMUT_R;
   float  b02 = dot(GAMUT_RB,v2);
   float  b12 = dot(GAMUT_RG,v2);
   bu = (GAMUT_B11*b02 - GAMUT_B01*b12) * GAMUT_BK;
   bv = (GAMUT_B00*b12 - GAMUT_B01*b02) * GAMUT_BK;
}

// test if we're inside gamut
int baryingamut(float bu,bv) {
   return (bu>=0 && bv>=0 && (bu+bv)<=1);
}
int xyingamut(vector xy) {
   float bu,bv; barysrgb(xy,bu,bv);
   return baryingamut(bu,bv);
}

// Gamut Mapping: "Absolute Colorimetric"
// Interpreted as projecting onto the closest edge of gamut triangle. 
// This implementation uses barycentric coords for the "sidedness" tests
// of the sample location (if out-of-gamut) before doing the projection.

// Version 1:
vector xytogamut_absolute_barycentric(vector xy) {
   float bu,bv; barysrgb(xy,bu,bv);
   if(baryingamut(bu,bv)) return xy;

   // this is fast because it uses the barycentric coords directly, but
   // not very accurate because the projections are not orthogonal to
   // the gamut's edges -- the result is that the chroma will be skewed 
   // more than it really should.
   float uvsum = bu+bv;
   if(uvsum>1) { bu/=uvsum; bv/=uvsum; }
   return GAMUT_R + GAMUT_RB*clamp(bu,0,1) + GAMUT_RG*clamp(bv,0,1);
}

// Version 2:
vector xytogamut_absolute_orthogonal(vector xy) {
   float bu,bv; barysrgb(xy,bu,bv);
   if(baryingamut(bu,bv)) return xy;

   // this still uses the bary coords as discriminants for which edge to
   // project onto, but then carries out an actual projection (unlike v1
   // which uses the non-orthogonal frame of the bary coords directly).
   if(bu<0) // red-green halfspace
      return GAMUT_R + GAMUT_NRG*clamp(dot(GAMUT_NRG,xy-GAMUT_R),0,GAMUT_LRG);

   if(bv<0) // red-blue halfspace
      return GAMUT_R + GAMUT_NRB*clamp(dot(GAMUT_NRB,xy-GAMUT_R),0,GAMUT_LRB);

   // green-blue halfspace
   return GAMUT_G + GAMUT_NGB*clamp(dot(GAMUT_NGB,xy-GAMUT_G),0,GAMUT_LGB);
}

// Pick one of the "Absolute Colorimetric" versions:
#define xytogamut_absolute   xytogamut_absolute_orthogonal

// Gamut Mapping: "Relative Colorimetric"
// Interpreted as moving toward the white point until the color is at
// the gamut triangle's edge.
// This implementation uses line-line intersections for the tests.
vector xytogamut_relative_nt(vector xy) {
   float  u;
   vector wp = WP_D65*{1,1,0};
   vector ws = xy - wp;
   // test against the red-green segment
   float det = (GAMUT_RG.y*ws.x) - (GAMUT_RG.x*ws.y);
   if(det>0) {
      u = ( (ws.x*GAMUT_RW.y) - (ws.y*GAMUT_RW.x) ) / det;
      if(u>=0 && u<=1) return GAMUT_R + GAMUT_RG*u;
   }
   // test against the red-blue segment
   det = (GAMUT_RB.x*ws.y) - (GAMUT_RB.y*ws.x);
   if(det>0) {
      u = ( (ws.y*GAMUT_RW.x) - (ws.x*GAMUT_RW.y) ) / det;
      if(u>=0 && u<=1) return GAMUT_R + GAMUT_RB*u;
   }
   // test against the green-blue segment
   det = (GAMUT_GB.y*ws.x) - (GAMUT_GB.x*ws.y);
   if(det>0) {
      u = ( (ws.x*GAMUT_GW.y) - (ws.y*GAMUT_GW.x) ) / det;
      if(u>=0 && u<=1) return GAMUT_G + GAMUT_GB*u;
   }
   return xy;
}
vector xytogamut_relative(vector xy) {
   if(xyingamut(xy)) return xy;
   return xytogamut_relative_nt(xy);
}


#define TOGAMUT_NONE       0

#define TOGAMUT_CLAMP      1
#define TOGAMUT_ABSOLUTE   2
#define TOGAMUT_RELATIVE   3

#define TOGAMUT_LAST       4


// Gamut-mapping dispachers for RGB colors
vector rgbtogamut(vector col; int methodid) {
   if(!methodid || methodid>=TOGAMUT_LAST) 
      return col;
   if(methodid==TOGAMUT_CLAMP) 
      return max(VZERO,col);

   vector c = col; float csign = 1;
   if(c.r<0 && c.g<0 && c.b<0) {
      csign = -1;
      c = abs(col);
   }
   vector xyy = rgbtoxyy(c);
   vector xy  = xyy * {1,1,0};
   vector ing,rgb;
   if(methodid==TOGAMUT_ABSOLUTE) ing = xytogamut_absolute(xy);
      else ing = xytogamut_relative(xy);
   ing = set(ing.x,ing.y,xyy.z);
   rgb = xyytorgb(ing);
   return rgb * csign * luma(col) / luma(rgb);
}

// XYZ -> In-Gamut Linear RGB
vector xyztorgbgamut(vector xyz) {
   vector rgb = xyztorgb(xyz);
   return rgbtogamut(rgb,TOGAMUT_RELATIVE);
}

// sRGB Gamma based on the sRGB specification, a copy of which can be found 
// at: www.color.org/sRGB.pdf, and the details in the original sRGB website
// can still be found via the waybackmachine at:
// http://web.archive.org/web/20030212204955/www.srgb.com/basicsofsrgb.htm

// sRGB -> Linear RGB
float fromSRGB(float srgb) {
   float c1 = 1.0 + 0.055;
   float c2 = 2.0 + 0.4;
   float out = srgb;
   if(srgb <= 0.04045) {
      out /= 12.92;
   } else {
      out = pow((out + 0.055) / c1, c2);
   }
   return out;
}
vector fromSRGB(vector srgb) {
   vector out;
   int i; for(i=0;i<3;i++) out[i] = fromSRGB(srgb[i]);
   return out;
}
vector4 fromSRGB(vector4 srgb) {
   vector c = fromSRGB((vector)srgb);
   return set(c.x,c.y,c.z,srgb.w);
}

// Linear RGB -> sRGB
float toSRGB(float rgb) {
   float c1 = 1.0 + 0.055;
   float c2 = 1.0 / 2.4;
   float out = rgb;
   if(rgb <= 0.0031308) {
      out *= 12.92;
   } else {
      out = c1 * pow(out,c2) - 0.055;
   }
   return out;
}
vector toSRGB(vector rgb) {
   vector out;
   int i; for(i=0;i<3;i++) out[i] = toSRGB(rgb[i]);
   return out;
}
vector4 toSRGB(vector4 rgb) {
   vector c = toSRGB((vector)rgb);
   return set(c.x,c.y,c.z,rgb.w);
}



//--------------------------------------------------------------------------------
// RGB
//--------------------------------------------------------------------------------

// Saturation
vector ccSat(vector col; float sat) {
   vector xyy = rgbtoxyy(col);
   vector xy  = xyy*{1,1,0};
   int wasok  = xyingamut(xy);
   vector wp  = WP_D65*{1,1,0};
   xy = wp + sat*(xy-wp);
   if(!xyingamut(xy)) {
      if(!wasok && sat>0) return col;
      xy = xytogamut_relative(xy);
   }
   vector rgb = xyytorgb(set(xy.x,xy.y,xyy.z));
   return rgb * luma(col) / luma(rgb);
}

// Hue Rotation Preserving Luminance
// Rotation is normalized, so [0,1] -> [0,360 degs]
vector ccHue(vector col; float rot) {
   vector c=col;
   if(rot!=0.) {
      vector l = LUMA_Default;
      
      //set x/y rotation matrices for grey vector
      float xrs=C_1_SQRT2;
      float xrc=xrs;
      matrix3 mat = set(1,0,0, 0,xrc,xrs, 0,-xrs,xrc);
      float yrs=-C_1_SQRT3;
      float yrc=C_SQRT2/C_SQRT3;
      mat *= set(yrc,0,-yrs, 0,1,0, yrs,0,yrc);
      
      //shear space to make the luminance plane horizontal
      vector ptmp=l*mat;  
      float dx=ptmp.x/ptmp.z;
      float dy=ptmp.y/ptmp.z;
      mat *= set(1,0,dx, 0,1,dy, 0,0,1);
      
      //rotate the hue
      float angle=rot*C_PI2;
      float zrs=sin(angle);
      float zrc=cos(angle);
      mat *= set(zrc,zrs,0, -zrs,zrc,0, 0,0,1);
      
      //unshear
      mat *= set(1,0,-dx, 0,1,-dy, 0,0,1);
      
      //un-rotate
      mat *= set(yrc,0,yrs, 0,1,0, -yrs,0,yrc);
      mat *= set(1,0,0, 0,xrc,-xrs, 0,xrs,xrc);
      
      c = c*mat;
   }

   return c;

}

// Returns the complementary rgb color to the one given (retaining luminance)
vector colortoopacity(vector c) {
   return rgbtogamut(ccSat(c,-1),TOGAMUT_RELATIVE);
}

// Brightness
float ccBright(float col; float b) {
   return b!=1.0 ? col*b : col;
}
vector ccBright(vector col; float b) {
   return b!=1.0 ? col*b : col;
}
vector ccBright(vector col; vector b) {
   return b!=VONE ? col*b : col;
}


// Contrast (pivot = 1)
float ccContrast(float val, cont) {
   return (val-1.0)*cont + 1.0;
}
vector ccContrast(vector val; float cont) {
   return (val-VONE)*cont + VONE;
}
vector ccContrast(vector val, cont) {
   return (val-VONE)*cont + VONE;
}


// Gamma
float ccGamma(float val, gam) {
   return pow(val,1.0/gam);
}
vector ccGamma(vector val; float gam) {
   return pow(val,1.0/gam);
}
vector ccGamma(vector val, gam) {
   return pow(val,VONE/gam);
}


// HSV Correction
vector ccHSV(vector col, hsv) {
   vector c = col;
   if(hsv.x!=0) c = ccHue(c,hsv.x);
   if(hsv.y!=1) c = ccSat(c,hsv.y);
   if(hsv.z!=1) c = ccBright(c,hsv.z);
   return c;
}



vector tonemap(vector rgb; float avgloglum, key, burn) {
   float rhk = key / avgloglum;
   float rhb = (burn!=0) ? pow(2.0, burn) : 1.0;
   float Lp = luma(rgb)*rhk;
   return rgb * Lp * (1.0+Lp*rhb) / (1.0+Lp);
}
vector tonemap(vector rgb; float avgloglum; vector key, burn) {
   vector rhk = key / avgloglum;
   vector rhb = set(
         (burn.x!=0) ? pow(2.0, burn.x) : 1.0,
         (burn.y!=0) ? pow(2.0, burn.y) : 1.0,
         (burn.z!=0) ? pow(2.0, burn.z) : 1.0 );
   vector Lp = luma(rgb)*rhk;
   return rgb * Lp * (1.0+Lp*rhb) / (1.0+Lp);
}


// Full CC
vector cc(
      vector col;
      vector c_hsv, c_contrast, c_bias, c_gain, c_gamma, c_tint;
      int    legalid; // 0=untouched, 1=clamp, 2=in-gamut
   ) 
{
   vector out = col;
   if(c_hsv!={0,1,1})      out = ccHSV(out,c_hsv);
   if(c_contrast!=VONE)    out = ccContrast(out,c_contrast);
   if(c_bias!=(vector)0.5) out = bias(out,c_bias);
   if(c_gain!=(vector)0.5) out = gain(out,c_gain);
   if(c_gamma!=VONE)       out = ccGamma(out,c_gamma);
   if(c_tint!=VONE)        out*= c_tint;
   if(legalid==1)          out = max(VZERO,out);
   return out;
}


// Full CC  -- This version with toggles is for the VOP wrapper
vector cctoggle (
      vector col;
      int    c_gamut;  // 0=none, 1=clamp, 2=absolute, 3=relative
      float  c_hrot; 
      vector c_bias, c_gain, c_gamma, c_cont;
      float  c_sat, c_val;
      vector c_tint;
      int    dohrot;
      int    dobias, dogain, dogamma, docont;
      int    dosat, doval;
      int    dotint;
   )
{
   vector out = col;
   vector tmp;
   if(dohrot && c_hrot!=0)
      out = ccHue(out,c_hrot);
   if(dobias && c_bias!=(vector)0.5) 
      out = bias(out,c_bias);
   if(dogain && c_gain!=(vector)0.5) 
      out = gain(out,c_gain);
   if(dogamma && c_gamma!=VONE)      
      out = ccGamma(out,c_gamma);
   if(docont && c_cont!=VONE)    
      out = ccContrast(out,c_cont);
   if(dosat && c_sat!=1)
      out = ccSat(out,c_sat);
   if(doval && c_val!=1)    
      out = ccBright(out,c_val);
   if(dotint && c_tint!=VONE)        
      out*= c_tint;
   return rgbtogamut(out,c_gamut);
}

#endif // End pyro_shade_h
