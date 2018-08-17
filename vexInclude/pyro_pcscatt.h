/****************************************************************************
          File:  pyro_pcscatt.h
       Created:  Jan 08, 2012
        Author:

   Description:  Support for PC-based scattered emission.
                 A little bit of baggage from pyro1 (updated and faster) that 
                 we'll keep around while a more robust renderer-side solution
                 is being worked on...

*****************************************************************************/



#ifndef pyro_pcscatt_h__GUARD
#define pyro_pcscatt_h__GUARD

#include <pyro_utils.h>

// The vop api for writing the clouds
// All clouds written and read in object space
void pyro_pcwrite(string firefile; vector emitcol) {
   vector po = ptransform("space:object",getblurP(0));
   if(firefile!="" && max(emitcol)>1e-4) {
      float  a  = C_PI_2 * length(wo_vspace(Du(P)+Dv(P)+Dw(P))); // = 2*pi*(r/4)
      pcwrite(     firefile
                  ,"interpolate",   1
                  ,"P",             po
                  ,"A",             a
                  ,"emit",          emitcol
              );
   }
}


// Given a pc handle representing fire intensities, calculate incident
// radiance from the fire volume.
// This "fast" version, assigns a constant density (given in "dens") to
// the entire intervening medium between the sample point and the pc points.
// It DOES NOT trace this interval to find actual density or occlusion.
// Therefore, it does not account for shadows from any opaque objects
// embedded in the fire/smoke volumes and the attenuation due to intervening
// density is only very roughly estimated -- result: a very fast "glow" from
// the fire, but not much else (accent on "fast":).
vector pyro_scatter_fast(int handle; bsdf Fsmoke; vector po,wo,dlocal,davg) {
   vector out = 0;
   float  d;
   vector pcP,pcE,wi,att;
   while(pciterate(handle)) {
      pcimport(handle,"P",pcP);
      pcimport(handle,"point.distance",d);
      pcimport(handle,"emit",pcE);
      wi  = normalize(pcP-po);
      att = eval_bsdf(Fsmoke,wo,wi) / (d*d + 1);
      float dl = max(dlocal);
      float da = max(davg);
      float dens = dl;
      if(dl<da) {
         float dh = da *.65;
         dens = (da-dh) * exp(dl-da) + dh;
      }
      att*= exp(-d*dens);
      out += pcE * att;
   }
   return out;
}

// Given a pc handle representing fire intensities, calculate incident
// radiance from the fire volume.
// This version traces for density/occlusion for each sample, so it's
// more accurate, but also slower.
vector pyro_scatter_trace(int handle; bsdf Fsmoke; vector po,wo) {
   vector out = 0;
   float  d;
   vector pcP,pcE,wi,hitOf,att;
   int    rv;
   while(pciterate(handle)) {
      pcimport(handle,"P",pcP);
      pcimport(handle,"point.distance",d);
      pcimport(handle,"emit",pcE);
      wi  = normalize(vtransform("space:object","space:camera", pcP - po ));
      att = eval_bsdf(Fsmoke,wo,wi) / (d*d + 1);
      rv  = trace(P, wi, 0,
               "raystyle" , "shadow" , 
               "maxdist"  , d        , 
               "Of"       , hitOf
            );
      if(rv>0) att *= VONE - clamp(hitOf,VZERO,VONE);
      out += pcE * att;
   }
   return out;
}


// Calc irradiance from a fire pc which was previously saved
// using pyro_pcwrite() above
vector pyro_scatter(string firefile, mode; 
                    bsdf Fsmoke; vector wo, dlocal,davg; float qual) 
{
   vector out = 0;
   vector po  = ptransform("space:object",getblurP(0));

   // calc the average area of the closest 10 pts
   int handle = pcopen(firefile,"P",{0,500,0},1e9,10);
   if(handle>=0) {
      float a, asum = 0;
      int npt = 0;
      while(pciterate(handle)) {
         pcimport(handle,"A",a);
         asum += a;
         npt++;
      }
      a = asum / max(1.0,(float)npt);  // ~ avg pt area
      npt = pcsize(handle);            // tot num pts in cloud
      pcclose(handle);

      // now open it using LOD based on avg pt-area and user-set "quality"
      int    fast   = (mode == "fast");
      float  q = qual; //fast ? qual : qual*0.5;
      float  thresh = a*2e4*exp(-5*q) + a;
      handle = pcopenlod(firefile,"P",po,8,
            "measure"       , "solidangle", 
            "area"          , "A"         , "threshold", thresh, 
            "aggregate:A"   , "sum"       , 
            "aggregate:emit", "sum"       , 
            "aggregate:P"   , "mean"
         );

      if(handle>=0) {
         out = fast ? pyro_scatter_fast(handle,Fsmoke,po,wo,dlocal,davg) :
                      pyro_scatter_trace(handle,Fsmoke,po,wo);
         pcclose(handle);
      }
   }

   return out*dPdz*.025*davg;

}
   
#endif // End pyro_pcscatt_h
