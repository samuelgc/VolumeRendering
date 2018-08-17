//
// This public-domain source is included in this distribution with the kind
// permission of its creator, Mario Marengo of AXYZ FX. For more information
// regarding this shader, please visit:
//
//	http://odforce.net/forum/index.php?showtopic=2262
//	http://www.sidefx.com/index.php?option=com_wrapper&Itemid=146
//

#ifndef __pcscatter__
#define __pcscatter__

#include "shading.h"

float vop_hgphase(float g; vector Wi, Wo) {
   float g2=g*g;
        return (1. - g2) / pow(1. + g2 - (2.*g*dot(-Wi,Wo)), 1.5);
}

float vop_hgphaseN(float g; vector Wi, Wo) {
   float g2=g*g;
   float num = g>0 ? -1.+g : 1.+g;
        return pow(num*num,1.5) / pow(1.+g2-(2.*g*dot(-Wi,Wo)),1.5);
}

float vop_ssBounceAtten(vector No,Ni,Li) {
   return 1.0 - ((1.0-dot(No,Ni))* (1.0-dot(No,Li)) / 2.0);
}

float vop_cdfSingle(float y) {
   float yy = clamp(y,0.0,1.0);
   return y<=0. ? 0. : (y>=1. ? 1. : yy*(2.0 + (yy-2.0)*yy*yy));
}

#if defined(VOP_SURFACE) || defined(VOP_DISPLACE) || defined(VOP_FOG)
void vop_pcIllum (int handle; string att) {
   vector p, n;
   vector illum;
   int status;
   while (pcunshaded(handle, att)) {
      pcimport(handle, "P", p); p = ow_space(p);
      pcimport(handle, "N", n); n = normalize(ow_nspace(n));
      illum = 0;
      illuminance(p, n, M_PI/2, LIGHT_DIFFUSE, "lightexport", "") {
         shadow(Cl);
         illum += Cl * diffuseBRDF(normalize(L), n);
      }
      status = pcexport(handle, att, illum);
   }
}
#endif

#if defined(VOP_SURFACE) || defined(VOP_DISPLACE) || defined(VOP_LIGHT) || \
    defined(VOP_FOG)
void vop_pcIrrad (int handle; string att) {
   vector p, n;
   vector irrad;
   int status;
   while (pcunshaded(handle, att)) {
      pcimport(handle, "P", p); p = ow_space(p);
      pcimport(handle, "N", n); n = normalize(ow_nspace(n));
      irrad = irradiance(p,n);
      status = pcexport(handle, att, irrad);
   }
}
#endif

float vop_icdfSingle(float x) {
   float rslt = 0.;
   if(x>=1.) {
      rslt = 1.;
   } else if(x>0.) {
      float a = pow(9. - (9.*x) + (1.73205080756887729353 *
                  sqrt(11.-(6.*x)-(21.*x*x)+(16.*x*x*x))), 0.333333333333333);
      float A = (2.*1.587401051968199*(1.-x)) / (1.442249570307408*a);
      float B = (1.259921049894873*a) / 2.080083823051904;
      float C = sqrt(1.0 + A + B);
      rslt = 0.5 * ( 1.0 + sqrt(2.0 - A - B + (2.0/C)) - C );
   }
   return rslt;
}


#if defined(VOP_SURFACE) || defined(VOP_DISPLACE) || defined(VOP_FOG)
vector vop_ssIntegMulti (
   string pcmap;
   vector Rdo;
   float sd;
   float bounce;
   int t_rgb;
   vector pcP;
   vector pcN;
   )
{
   vector Xi,Ni;
   vector Xo = pcP;
   vector No = normalize(pcN);
   vector ld = Rdo*sd;
   float ld1 = max(ld);
   int handle = pcopen(pcmap, "P", Xo, ld1, (int)1e9);
   vop_pcIllum(handle,"illum");
   float r,ptarea;
   vector ssm=0, ptillum=0;
   while (pciterate(handle)) {
      pcimport(handle, "P", Xi);
      pcimport(handle, "N", Ni);
      pcimport(handle, "point.distance", r);
      pcimport(handle, "ptarea", ptarea);
      pcimport(handle, "illum", ptillum);
      Ni = normalize(Ni);
      vector Li = (Xo-Xi)/ld1;
      float kb = vop_ssBounceAtten(No,Ni,Li);
      kb = lerp(1.0,kb,bounce);
      if(kb>0.0 ) {
         if(t_rgb)
         {
            int wave;
            for(wave=0;wave<3;wave++) {
               setcomp( ssm,
                        getcomp(ssm,wave) +
                           kb * getcomp(ptillum,wave) * ptarea *
                           (1-smooth(0,getcomp(ld,wave),r)),
                        wave
                      );
            }
         }
         else
            ssm += kb * ptillum * ptarea * (1-smooth(0,ld1,r));
      }
   }
   pcclose(handle);
   if(!t_rgb) ssm*=Rdo;
   float norm = 3.0*ld1*ld1*M_PI / 10.0;
   return ssm / norm;
}


vector vop_ssIntegSingle (
   vector Rd;
   float sd;
   float g;
   float eta;
   int samples;
   float tbias;
   int t_rgb;
   vector PP;
   vector NN;
   vector II;
   )
{
   float Kro, Kto, Kri, Kti;
   vector Xi, Wi, Ni;
   int samp;
   vector Psamp;
   float phase,WiNi,AWiNi,spi,spo,ksss;
   string oname = getobjectname();
   vector Xo = PP;
   vector No = normalize(NN);
   vector Wo = -normalize(II);
   vector lu = Rd*sd;
   float lu1 = max(lu);
   float ieta = 1.0 / eta;
   float ieta2 = ieta * ieta;
   vector To = normalize(refract(-Wo,No,ieta));
   vector Wpo = -To;
   float gg = clamp(g,-0.998,.998);
   if(eta!=1.0) { fresnel(Wpo,-No,ieta,Kro,Kto); }
      else Kto = 1.0;
   vector scatt = 0;
   vector realsamples = 0;
   float hitD = rayhittest(Xo,To*1e6,tbias,"scope",oname);
   float spoMax = hitD<0. ? lu1 : min(hitD,lu1);
   float terr = tbias;
   float sinc = (1.0-2.0*terr)/(float)(samples);
   float ss = terr;
   float ssbase= ss;
   if(t_rgb) {
      vector maxadj = set(
            vop_cdfSingle(clamp(spoMax/lu.x,0.,1.)),
            vop_cdfSingle(clamp(spoMax/lu.y,0.,1.)),
            vop_cdfSingle(clamp(spoMax/lu.z,0.,1.))
         );
      int wave;
      for(wave=0;wave<3;wave++)
      {
         ss=ssbase=terr;
         float luk = getcomp(lu,wave);
         for(samp=0; samp<samples; samp++)
         {
            ss = ssbase+sinc*nrandom();
            spo = spoMax*vop_icdfSingle(ss*getcomp(maxadj,wave));
            ssbase+=sinc;
            Psamp = Xo + (To * spo);
            illuminance(Psamp, No, M_PI, LIGHT_DIFFUSE, "lightexport", "")
            {
               Wi = normalize(L);
               hitD = rayhittest(Psamp,L,Xi,Ni,0.,"scope",oname);
               if(hitD>0.) {
                  setcomp(realsamples,getcomp(realsamples,wave)+1,wave);
                  Ni = normalize(Ni);
                  WiNi = dot(Wi,Ni);
                  AWiNi = abs(WiNi);
                  spi = distance(Psamp,Xi) * AWiNi /
                              sqrt(1.0 - ieta2 * (1.0 - AWiNi*AWiNi));
		  if(spi <= luk && WiNi>0.) {
		      if(eta!=1.0) { fresnel(-Wi,Ni,ieta,Kri,Kti); }
		      else Kti = 1.0;
		      float f = Kti * Kto;
		      phase = vop_hgphaseN(gg,Wi,Wpo);
		      ksss = f * phase * (1-smooth(0,luk,spi));
		      setcomp(scatt,getcomp(scatt,wave) +
			      getcomp(Cl,wave)*WiNi*ksss,
			      wave);
                  }
               }
            }
         }
      }
   } else {
      ss=ssbase=terr;
      float icdfmax = vop_cdfSingle(clamp(spoMax/lu1,0.,1.));
      for(samp=0; samp<samples; samp++)
      {
         ss = ssbase+nrandom()*sinc;
         spo = spoMax*vop_icdfSingle(ss*icdfmax);
         ssbase += sinc;
         Psamp = Xo + (To * spo);
         illuminance(Psamp, No, M_PI, LIGHT_DIFFUSE, "lightexport", "")
         {
            Wi = normalize(L);
            hitD = rayhittest(Psamp,L,Xi,Ni,0.,"scope",oname);
            if(hitD>0.) {
               realsamples += 1;
               Ni = normalize(Ni);
               WiNi = dot(Wi,Ni);
               AWiNi = abs(WiNi);
               spi = distance(Psamp,Xi) * AWiNi /
                           sqrt(1.0 - ieta2 * (1.0 - AWiNi*AWiNi));
               if(spi <= lu1 && WiNi>0.) {
                  if(eta!=1.0) { fresnel(-Wi,Ni,ieta,Kri,Kti); }
                     else Kti = 1.0;
                  float f = Kti * Kto;
                     phase = vop_hgphaseN(gg,Wi,Wpo);
                     ksss = f * phase * (1-smooth(0,lu1,spi));
                  scatt+= Cl * Rd * WiNi * ksss;
               }
            }
         }
      }
   }
      return 2.0 * scatt / realsamples;
}
#endif

#endif // End sss_GUARD__
