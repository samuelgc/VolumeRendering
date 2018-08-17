/****************************************************************************
   Description:  Noise functions for pyro.
                 They work in all contexts. 
                 They extend the noises in vex/include/voplib.h

                  Input: Pos 1 to 4 dimensions
                 Output: float or vector in [0,1]

                 -Normalized over all dimensions, octaves, and roughness
                 -Continuous octaves
                 -Frequency (lacunarity) and amplitude (roughness)
                 -Space warping (variable lacunarity) and gradient
                  warping now part of the basis definition, not
                  the fractal construction -- much cleaner and now
                  available for all noise and fractal types as a result.
                 -Added support for simplex noise (in all its native variants)
                  as well as its built-in gradient (so local vex-based
                  gradient calculation is skipped for the xnoise() basis and
                  the built-in xnoised() function used instead for speed).
                  NOTE: None of the xnoise() functions are documented!
                 -All noises and fractals have a unified API.
                 -No assumptions made about the filter size. Left to the
                  caller. Cleaner and context-agnostic.
                 -All noise bases support all dimensions. This is forced so the
                  API is consistent. If the underlying native implementation
                  of a basis doesn't support a certain dimension, then the
                  wrapper here will simply pad it or build it up to at least
                  not fail. Therefore, all position, frequency, offset, and 
                  all space-related quantities are 4-dimensional args.

*****************************************************************************/

#ifndef pyro_noise2_h__GUARD
#define pyro_noise2_h__GUARD
   

#include <pyro_utils.h>


//-------------------------------------------------------------------------------
// Noise Statistics
// Simple stats (min, max, average) on the available VEX noises.
// Uses a variance threshold (applied to the median of the running min 
// and max values) to determine when to stop generating samples.

struct nsdata {
   string name;
   float  min,max,avg;
   int    symmetric;
}

float stats_samplef(string name; vector4 p; int dim; float period) 
{
   int per = (int)period;
   int seed; float f1,f2;
   if(name=="perlin") {
      if(dim==1) return float(noise(p.x));
      else if(dim==2) return float(noise(p.x,p.y));
      else if(dim==3) return float(noise((vector)p));
      else return float(noise(p));
   } else if(name=="pperlin") {
      if(dim==1) return float(pnoise(p.x,per));
      else if(dim==2) return float(pnoise(p.x,p.y,per,per));
      else if(dim==3) return float(pnoise((vector)p,(vector)per));
      else return float(pnoise(p,(vector4)per));
   } else if(name=="simplex") {
      if(dim==1) return float(xnoise(p.x));
      else if(dim==2) return float(xnoise(p.x,p.y));
      else if(dim==3) return float(xnoise((vector)p));
      else return float(xnoise(p));
   } else if(name=="sparse") {
      if(dim==1) return float(snoise(set(p.x,0,0)));
      else if(dim==2) return float(snoise(set(p.x,p.y,0)));
      else if(dim==3) return float(snoise((vector)p));
      else return float(snoise((vector)p+set(p.w,-p.w,avg((vector)p))));
   } else if(name=="flow") {
      if(dim==1) return float(flownoise((vector)p.x,0));
      else if(dim==2) return float(flownoise(p.x,p.y,0));
      else if(dim==3) return float(flownoise((vector)p,0));
      else return float(flownoise(p,0));
   } else if(name=="pflow") {
      if(dim==1) return float(flowpnoise((vector)p.x,per,0));
      else if(dim==2) return float(flowpnoise(p.x,p.y,per,per,0));
      else if(dim==3) return float(flowpnoise((vector)p,(vector)per,0));
      else return float(flowpnoise(p,(vector4)per,0));
   } else if(name=="worleyFA") {
      if(dim==1) {
         wnoise(p.x,seed,f1,f2); 
         return f1;
      } else if(dim==2) {
         wnoise(p.x,p.y,seed,f1,f2); 
         return f1;
      } else if(dim==3) {
         wnoise((vector)p,seed,f1,f2);
         return f1;
      } else {
         wnoise(p,seed,f1,f2); 
         return f1;
      }
   } else if(name=="worleyFB") {
      if(dim==1) {
         wnoise(p.x,seed,f1,f2); 
         return (f2-f1);
      } else if(dim==2) {
         wnoise(p.x,p.y,seed,f1,f2); 
         return (f2-f1);
      } else if(dim==3) {
         wnoise((vector)p,seed,f1,f2);
         return (f2-f1);
      } else {
         wnoise(p,seed,f1,f2);
         return (f2-f1);
      }
   } else if(name=="mworleyFA") {
      if(dim==1) {
         mwnoise(p.x,seed,f1,f2); 
         return f1;
      } else if(dim==2) {
         mwnoise(p.x,p.y,seed,f1,f2); 
         return f1;
      } else if(dim==3) {
         mwnoise((vector)p,seed,f1,f2);
         return f1;
      } else {
         mwnoise(p,seed,f1,f2); 
         return f1;
      }
   } else if(name=="mworleyFB") {
      if(dim==1) {
         mwnoise(p.x,seed,f1,f2); 
         return (f2-f1);
      } else if(dim==2) {
         mwnoise(p.x,p.y,seed,f1,f2); 
         return (f2-f1);
      } else if(dim==3) {
         mwnoise((vector)p,seed,f1,f2);
         return (f2-f1);
      } else {
         mwnoise(p,seed,f1,f2);
         return (f2-f1);
      }
   } else if(name=="cworleyFA") {
      if(dim==1) {
         cwnoise(p.x,seed,f1,f2); 
         return f1;
      } else if(dim==2) {
         cwnoise(p.x,p.y,seed,f1,f2); 
         return f1;
      } else if(dim==3) {
         cwnoise((vector)p,seed,f1,f2);
         return f1;
      } else {
         cwnoise(p,seed,f1,f2); 
         return f1;
      }
   } else if(name=="cworleyFB") {
      if(dim==1) {
         cwnoise(p.x,seed,f1,f2); 
         return (f2-f1);
      } else if(dim==2) {
         cwnoise(p.x,p.y,seed,f1,f2); 
         return (f2-f1);
      } else if(dim==3) {
         cwnoise((vector)p,seed,f1,f2);
         return (f2-f1);
      } else {
         cwnoise(p,seed,f1,f2);
         return (f2-f1);
      }
   } else if(name=="alligator") {
      if(dim==1) return float(anoise(set(p.x,0,0)));
      else if(dim==2) return float(anoise(set(p.x,p.y,0)));
      else if(dim==3) return float(anoise((vector)p));
      else return float(anoise((vector)p+set(p.w,-p.w,avg((vector)p))));
   }

   return 0;
}


vector stats_samplev(string name; vector4 p; int dim; float period) 
{
   int per = (int)period;
   int seed; float f1,f2;
   if(name=="perlin") {
      if(dim==1) return vector(noise(p.x));
      else if(dim==2) return vector(noise(p.x,p.y));
      else if(dim==3) return vector(noise((vector)p));
      else return vector(noise(p));
   } else if(name=="pperlin") {
      if(dim==1) return vector(pnoise(p.x,per));
      else if(dim==2) return vector(pnoise(p.x,p.y,per,per));
      else if(dim==3) return vector(pnoise((vector)p,(vector)per));
      else return vector(pnoise(p,(vector4)per));
   } else if(name=="simplex") {
      if(dim==1) return vector(xnoise(p.x));
      else if(dim==2) return vector(xnoise(p.x,p.y));
      else if(dim==3) return vector(xnoise((vector)p));
      else return vector(xnoise(p));
   } else if(name=="sparse") {
      if(dim==1) return vector(snoise(set(p.x,0,0)));
      else if(dim==2) return vector(snoise(set(p.x,p.y,0)));
      else if(dim==3) return vector(snoise((vector)p));
      else return vector(snoise((vector)p+set(p.w,-p.w,avg((vector)p))));
   } else if(name=="flow") {
      if(dim==1) return vector(flownoise((vector)p.x,0));
      else if(dim==2) return vector(flownoise(p.x,p.y,0));
      else if(dim==3) return vector(flownoise((vector)p,0));
      else return vector(flownoise(p,0));
   } else if(name=="pflow") {
      if(dim==1) return vector(flowpnoise((vector)p.x,per,0));
      else if(dim==2) return vector(flowpnoise(p.x,p.y,per,per,0));
      else if(dim==3) return vector(flowpnoise((vector)p,(vector)per,0));
      else return vector(flowpnoise(p,(vector4)per,0));
   } else if(name=="worleyFA") {
      return 0;
   } else if(name=="worleyFB") {
      return 0;
   } else if(name=="mworleyFA") {
      return 0;
   } else if(name=="mworleyFB") {
      return 0;
   } else if(name=="cworleyFA") {
      return 0;
   } else if(name=="cworleyFB") {
      return 0;
   } else if(name=="alligator") {
      if(dim==1) return vector(anoise(set(p.x,0,0)));
      else if(dim==2) return vector(anoise(set(p.x,p.y,0)));
      else if(dim==3) return vector(anoise((vector)p));
      else return vector(anoise((vector)p+set(p.w,-p.w,avg((vector)p))));
   }

   return 0;
}


void nstats(float varthresh, spacesize,freq; int minsamps, maxsamps, repeats) 
{
   string names[] = {
      "perlin", "pperlin", "simplex", "sparse", "flow", "pflow",
      "worleyFA", "worleyFB", "mworleyFA", "mworleyFB",
      "cworleyFA", "cworleyFB", "alligator" };

   string  rtnames[] = {"float","vector"};
   string  rtprefix[] = {"f","v"};
   float   fmin,fmax,fsum,favg;
   vector  vmin,vmax,vsum,vavg;
   vector4 p;

   float thresh2 = varthresh*varthresh;

   string name,rtname;
   int i,dim,rtype; 
   float ksize[] = array(1,0.5,1.0/3.0,1.0/4.0);
   float ssize;
   for(i=0;i<len(names);i++) 
   {
      for(dim=4;dim<=4;dim++) 
      {
         ssize = spacesize * ksize[dim-1];
         for(rtype=0;rtype<2;rtype++) 
         {
            name   = names[i];
            rtname = rtnames[rtype];
            //printf("// Processing: \"%s\", Dim: %g, Type: \"%s\"",name,dim,rtname);

            // sample vars
            float  fsamp, d = 1e9, sx, sy, sz;
            vector vsamp;
            fmin = 1e9; fmax = -1e9; favg = 0; fsum = 0;
            vmin = 1e9; vmax = -1e9; vavg = 0; vsum = 0;

            // variance vars
            float dprev  = 1e9, dmean = 0, dvsum = 0, dvar, dd;
            int   done = 0, count = 0, j, jp1, under;

            // variance loop (var based on distance between curr and prev samples)
            for(j=0;!done;j++) 
            {
               jp1 = j+1;

               // draw a random 4D position 'p' from a space of size 'spacesize'
               p = ( vector4(nrandom("twister")) - 0.5 ) * ssize * freq;

               // sample the current noise at 'p'
               // The distance 'd' between the running min and max values is what
               // we're trying to reduce the variance of.
               if(rtype) {
                  vsamp = stats_samplev(name,p,dim,ssize);
                  vsum += vsamp;
                  vavg = vsum / jp1;
                  vmin = min(vsamp,vmin);
                  vmax = max(vsamp,vmax);
                  d = 0.5 * (max(vmax)+min(vmin));
               } else {
                  fsamp = stats_samplef(name,p,dim,ssize);
                  fsum += fsamp;
                  favg = fsum / jp1;
                  fmin = min(fsamp,fmin);
                  fmax = max(fsamp,fmax);
                  d = 0.5 * (fmax+fmin);
               }

               // calculate the running mean 'dmean' and variance 'dvar' of the
               // distance 'd' between running min and max
               dd     = d - dmean;
               dmean += dd / jp1;
               dvsum += dd * (d - dmean);
               dvar   = dvsum / j;

               // check if under threshold and stop if we are
               under = (dvar < thresh2);
               if(under) ++count; else count = 0;
               if(jp1<maxsamps) 
                  done = jp1>=minsamps && count>repeats;
               else
                  done = jp1>=minsamps;
            }

            // Report results
            //printf(",\t Samples: %d\n",j);
            float hvar = sqrt(dvar) * 0.5;
            float lo,hi,av;
            if(rtype) {
               lo = min(vmin); if(lo>0) lo = max(0,lo-hvar);
               hi = max(vmax)+hvar;
               av = efit(avg(vavg),lo,hi,0,1);
               printf("#define ns_%s%s%d  nsdata(%g,%g,%g) // +/- %g [%d]\n",
                     rtprefix[rtype],name,dim, lo,hi,av,hvar,j);
            } else {
               lo = fmin-hvar;
               hi = fmax+hvar;
               av = efit(favg,lo,hi,0,1);
               printf("#define ns_%s%s%d  nsdata(%g,%g,%g) // +/- %g [%d]\n",
                     rtprefix[rtype],name,dim, lo,hi,av,hvar,j);
            }
         }
      }
   }
}



//-------------------------------------------------------------------------------
// The constants used in this module were generated with the call:
//
//    nstats(0.005, 10,20, 1000, (int)5e6, 0);
//
// which produces the following output:

#define ns_fperlin1    \
   nsdata ( "perlin"    , 0.248834   , 0.767147 , 0.488376  , 1 ) // +/- 0.0024
#define ns_vperlin1    \
   nsdata ( "perlin"    , 0.229499   , 0.761877 , 0.509742  , 1 ) // +/- 0.0024
#define ns_fperlin2    \
   nsdata ( "perlin"    , 0.136616   , 0.864593 , 0.499529  , 1 ) // +/- 0.0024
#define ns_vperlin2    \
   nsdata ( "perlin"    , 0.033912   , 0.946057 , 0.511009  , 1 ) // +/- 0.0025
#define ns_fperlin3    \
   nsdata ( "perlin"    , 0.0832587  , 0.920337 , 0.497705  , 1 ) // +/- 0.0024
#define ns_vperlin3    \
   nsdata ( "perlin"    , 0.0013614  , 0.996832 , 0.335704  , 1 ) // +/- 0.0032
#define ns_fpperlin1   \
   nsdata ( "pperlin"   , 0.272539   , 0.732348 , 0.493888  , 1 ) // +/- 0.0024
#define ns_vpperlin1   \
   nsdata ( "pperlin"   , 0.24574    , 0.7608   , 0.493501  , 1 ) // +/- 0.0024
#define ns_fpperlin2   \
   nsdata ( "pperlin"   , 0.128858   , 0.839149 , 0.521924  , 1 ) // +/- 0.0024
#define ns_vpperlin2   \
   nsdata ( "pperlin"   , 0.0974457  , 0.914321 , 0.492854  , 1 ) // +/- 0.0024
#define ns_fpperlin3   \
   nsdata ( "pperlin"   , 0.0777629  , 0.911734 , 0.50605   , 1 ) // +/- 0.0024
#define ns_vpperlin3   \
   nsdata ( "pperlin"   , 0.0191109  , 0.982488 , 0.402437  , 1 ) // +/- 0.0025
#define ns_fsimplex1   \
   nsdata ( "simplex"   , 0.0135945  , 0.980643 , 0.503308  , 1 ) // +/- 0.0024
#define ns_vsimplex1   \
   nsdata ( "simplex"   , 0.00470505 , 0.979253 , 0.508709  , 1 ) // +/- 0.0024
#define ns_fsimplex2   \
   nsdata ( "simplex"   , 0.100222   , 0.909426 , 0.494677  , 1 ) // +/- 0.0024
#define ns_vsimplex2   \
   nsdata ( "simplex"   , 0.0576417  , 0.958953 , 0.491143  , 1 ) // +/- 0.0024
#define ns_fsimplex3   \
   nsdata ( "simplex"   , 0.15302    , 0.850784 , 0.497038  , 1 ) // +/- 0.0024
#define ns_vsimplex3   \
   nsdata ( "simplex"   , 0.0434933  , 0.970121 , 0.315176  , 1 ) // +/- 0.0076
#define ns_fsparse1    \
   nsdata ( "sparse"    , -1.05121   , 1.41258  , 0.437552  , 1 ) // +/- 0.0024
#define ns_vsparse1    \
   nsdata ( "sparse"    , -1.84633   , 1.41258  , 0.563931  , 1 ) // +/- 0.0025
#define ns_fsparse2    \
   nsdata ( "sparse"    , -1.85569   , 1.8013   , 0.514998  , 1 ) // +/- 0.0027
#define ns_vsparse2    \
   nsdata ( "sparse"    , -2.28436   , 2.08765  , 0.520132  , 1 ) // +/- 0.0025
#define ns_fsparse3    \
   nsdata ( "sparse"    , -2.34351   , 2.43843  , 0.49609   , 1 ) // +/- 0.0174
#define ns_vsparse3    \
   nsdata ( "sparse"    , -2.71525   , 2.64793  , 0.504632  , 1 ) // +/- 0.0256
#define ns_fflow1      \
   nsdata ( "flow"      , 0.191786   , 0.838335 , 0.476434  , 1 ) // +/- 0.0024
#define ns_vflow1      \
   nsdata ( "flow"      , 0.156953   , 0.847188 , 0.498596  , 1 ) // +/- 0.0024
#define ns_fflow2      \
   nsdata ( "flow"      , 0.110068   , 0.907473 , 0.489255  , 1 ) // +/- 0.0024
#define ns_vflow2      \
   nsdata ( "flow"      , 0.019407   , 0.977013 , 0.330133  , 1 ) // +/- 0.0032
#define ns_fflow3      \
   nsdata ( "flow"      , 0.0972697  , 0.879663 , 0.514725  , 1 ) // +/- 0.0024
#define ns_vflow3      \
   nsdata ( "flow"      , 0.0789278  , 0.909136 , 0.506907  , 1 ) // +/- 0.0024
#define ns_fpflow1     \
   nsdata ( "pflow"     , 0.192796   , 0.835272 , 0.483768  , 1 ) // +/- 0.0024
#define ns_vpflow1     \
   nsdata ( "pflow"     , 0.192727   , 0.834885 , 0.481012  , 1 ) // +/- 0.0024
#define ns_fpflow2     \
   nsdata ( "pflow"     , 0.0875699  , 0.872108 , 0.526021  , 1 ) // +/- 0.0025
#define ns_vpflow2     \
   nsdata ( "pflow"     , 0.0681927  , 0.928206 , 0.502054  , 1 ) // +/- 0.0025
#define ns_fpflow3     \
   nsdata ( "pflow"     , 0.0931273  , 0.896028 , 0.506575  , 1 ) // +/- 0.0024
#define ns_vpflow3     \
   nsdata ( "pflow"     , 0.0427369  , 0.940558 , 0.509313  , 1 ) // +/- 0.0024
#define ns_fworleyFA1  \
   nsdata ( "worley"    , 0          , 0.742495 , 0.0740117 , 0 ) // +/- 0.0024
#define ns_vworleyFA1  \
   nsdata ( "worley"    , 0          , 0.742495 , 0.0740117 , 0 ) // +/- 0.0024
#define ns_fworleyFA2  \
   nsdata ( "worley"    , 0          , 1.15271  , 0.108373  , 0 ) // +/- 0.0256
#define ns_vworleyFA2  \
   nsdata ( "worley"    , 0          , 1.15271  , 0.108373  , 0 ) // +/- 0.0256
#define ns_fworleyFA3  \
   nsdata ( "worley"    , 0          , 1.18895  , 0.159684  , 0 ) // +/- 0.0256
#define ns_vworleyFA3  \
   nsdata ( "worley"    , 0          , 1.18895  , 0.159684  , 0 ) // +/- 0.0256
#define ns_fworleyFB1  \
   nsdata ( "worley"    , 0          , 0.902963 , 0.118548  , 0 ) // +/- 0.0025
#define ns_vworleyFB1  \
   nsdata ( "worley"    , 0          , 0.902963 , 0.118548  , 0 ) // +/- 0.0025
#define ns_fworleyFB2  \
   nsdata ( "worley"    , 0          , 1.24931  , 0.108399  , 0 ) // +/- 0.0256
#define ns_vworleyFB2  \
   nsdata ( "worley"    , 0          , 1.24931  , 0.108399  , 0 ) // +/- 0.0256
#define ns_fworleyFB3  \
   nsdata ( "worley"    , 0          , 1.1101   , 0.118099  , 0 ) // +/- 0.0181
#define ns_vworleyFB3  \
   nsdata ( "worley"    , 0          , 1.1101   , 0.118099  , 0 ) // +/- 0.0181
#define ns_fmworleyFA1  \
   nsdata( "mworley"    , 0          , 0.587001 , 0.0971886 , 0 ) // +/- 0.0025
#define ns_vmworleyFA1  \
   nsdata( "mworley"    , 0          , 0.587001 , 0.0971886 , 0 ) // +/- 0.0025
#define ns_fmworleyFA2  \
   nsdata( "mworley"    , 0          , 1.29428  , 0.314845  , 0 ) // +/- 0.0059487
#define ns_vmworleyFA2  \
   nsdata( "mworley"    , 0          , 1.29428  , 0.314845  , 0 ) // +/- 0.0059487
#define ns_fmworleyFA3  \
   nsdata( "mworley"    , 0          , 1.56603    ,0.398481  , 0 ) // +/- 0.0124397
#define ns_vmworleyFA3  \
   nsdata( "mworley"    , 0          , 1.56603    ,0.398481  , 0 ) // +/- 0.0124397
#define ns_fmworleyFB1  \
   nsdata( "mworley"    , 0          , 0.618887 , 0.178215   , 0 ) // +/- 0.00249998
#define ns_vmworleyFB1  \
   nsdata( "mworley"    , 0          , 0.618887 , 0.178215   , 0 ) // +/- 0.00249998
#define ns_fmworleyFB2  \
   nsdata( "mworley"    , 0          , 1.25947  , 0.183175   , 0 ) // +/- 0.0153326
#define ns_vmworleyFB2  \
   nsdata( "mworley"    , 0          , 1.25947  , 0.183175   , 0 ) // +/- 0.0153326
#define ns_fmworleyFB3  \
   nsdata( "mworley"    , 0          , 1.45466    , 0.161667 , 0 ) // +/- 0.0213142
#define ns_vmworleyFB3  \
   nsdata( "mworley"    , 0          , 1.45466    , 0.161667 , 0 ) // +/- 0.0213142
#define ns_fcworleyFA1  \
   nsdata( "cworley"    , 0          , 0.587747   , 0.0978262, 0 ) // +/- 0.0025
#define ns_vcworleyFA1  \
   nsdata( "cworley"    , 0          , 0.587747   , 0.0978262, 0 ) // +/- 0.0025
#define ns_fcworleyFA2  \
   nsdata( "cworley"    , 0          , 0.901443   , 0.320862 , 0 ) // +/- 0.00593521
#define ns_vcworleyFA2  \
   nsdata( "cworley"    , 0          , 0.901443   , 0.320862 , 0 ) // +/- 0.00593521
#define ns_fcworleyFA3  \
   nsdata( "cworley"    , 0          , 0.843453   , 0.406956, 0 ) // +/- 0.00592824
#define ns_vcworleyFA3  \
   nsdata( "cworley"    , 0          , 0.843453   , 0.406956, 0 ) // +/- 0.00592824
#define ns_fcworleyFB1  \
   nsdata( "cworley"    , 0          , 0.619793   , 0.178844, 0 ) // +/- 0.0025
#define ns_vcworleyFB1  \
   nsdata( "cworley"    , 0          , 0.619793   , 0.178844, 0 ) // +/- 0.0025
#define ns_fcworleyFB2  \
   nsdata( "cworley"    , 0          , 0.849381   , 0.185333, 0 ) // +/- 0.00514269
#define ns_vcworleyFB2  \
   nsdata( "cworley"    , 0          , 0.849381   , 0.185333, 0 ) // +/- 0.00514269
#define ns_fcworleyFB3  \
   nsdata( "cworley"    , -0.0114286 , 0.776565   , 0.165661, 0 ) // +/- 0.0114286
#define ns_vcworleyFB3  \
   nsdata( "cworley"    , -0.0114286 , 0.776565   , 0.165661, 0 ) // +/- 0.0114286
#define ns_falligator1 \
   nsdata ( "alligator" , 0          , 0.897279 , 0.13911   , 0 ) // +/- 0.0024
#define ns_valligator1 \
   nsdata ( "alligator" , 0          , 0.931199 , 0.132454  , 0 ) // +/- 0.0024
#define ns_falligator2 \
   nsdata ( "alligator" , 0          , 0.981734 , 0.117792  , 0 ) // +/- 0.0025
#define ns_valligator2 \
   nsdata ( "alligator" , 0          , 0.980294 , 0.126717  , 0 ) // +/- 0.0024
#define ns_falligator3 \
   nsdata ( "alligator" , 0          , 0.993732 , 0.117951  , 0 ) // +/- 0.0032
#define ns_valligator3 \
   nsdata ( "alligator" , 0          , 0.992102 , 0.128566  , 0 ) // +/- 0.0025
#define ns_fperlin4 \
   nsdata ( "perlin"    , 0.0168713 , 0.998413  , 0.507642  , 1 ) // +/- 0.0073
#define ns_vperlin4 \
   nsdata ( "perlin"    , 0.00576016 , 1.025    , 0.518260  , 1 ) // +/- 0.0037
#define ns_fpperlin4   \
   nsdata ( "pperlin"   , 0.154528   , 0.828153 , 0.511577  , 1 ) // +/- 0.0024
#define ns_vpperlin4   \
   nsdata ( "pperlin"   , 0.149949   , 0.853128 , 0.49744   , 1 ) // +/- 0.0024
#define ns_fsimplex4   \
   nsdata ( "simplex"   , 0.0943673  , 0.912882 , 0.503625  , 1 ) // +/- 0.0064
#define ns_vsimplex4   \
   nsdata ( "simplex"   , 0.13602    , 0.848679 , 0.510355  , 1 ) // +/- 0.0025
#define ns_fsparse4    \
   nsdata ( "sparse"    , -2.18691   , 2.46426  , 0.476393  , 1 ) // +/- 0.0064
#define ns_vsparse4    \
   nsdata ( "sparse"    , -2.59173   , 2.50891  , 0.506553  , 1 ) // +/- 0.0145
#define ns_fflow4      \
   nsdata ( "flow"      , 0.0541632  , 0.942907 , 0.501736  , 1 ) // +/- 0.0025
#define ns_vflow4      \
   nsdata ( "flow"      , 0.0834745  , 0.893131 , 0.514653  , 1 ) // +/- 0.0024
#define ns_fpflow4     \
   nsdata ( "pflow"     , 0.144938   , 0.852499 , 0.501408  , 1 ) // +/- 0.0024
#define ns_vpflow4     \
   nsdata ( "pflow"     , 0.155242   , 0.840548 , 0.5022    , 1 ) // +/- 0.0024
#define ns_fworleyFA4  \
   nsdata ( "worley"    , 0          , 1.19425  , 0.314428  , 0 ) // +/- 0.0024
#define ns_vworleyFA4  \
   nsdata ( "worley"    , 0          , 1.19425  , 0.314428  , 0 ) // +/- 0.0024
#define ns_fworleyFB4  \
   nsdata ( "worley"    , 0          , 1.53913  , 0.1402    , 0 ) // +/- 0.0512
#define ns_vworleyFB4  \
   nsdata ( "worley"    , 0          , 1.53913  , 0.1402    , 0 ) // +/- 0.0512
#define ns_fmworleyFA4  \
   nsdata ( "mworley"   , 0.00495732 , 1.7116   , 0.482286  , 0 ) // +/- 0.0068835
#define ns_vmworleyFA4  \
   nsdata ( "mworley"   , 0.00495732 , 1.7116   , 0.482286  , 0 ) // +/- 0.0068835
#define ns_fmworleyFB4  \
   nsdata ( "mworley"   ,-0.0163645  , 1.42481  , 0.159796  , 0 ) // +/- 0.0163645 
#define ns_vmworleyFB4  \
   nsdata ( "mworley"   ,-0.0163645  , 1.42481  , 0.159796  , 0 ) // +/- 0.0163645
#define ns_fcworleyFA4  \
   nsdata ( "cworley"   , 0.0274073  , 0.690122 , 0.521913  , 0 ) // +/- 0.00249996
#define ns_vcworleyFA4  \
   nsdata ( "cworley"   , 0.0274073  , 0.690122 , 0.521913  , 0 ) // +/- 0.00249996
#define ns_fcworleyFB4  \
   nsdata ( "cworley"   , 0          , 0.647968 , 0.156623  , 0 ) // +/- 0.00587861
#define ns_vcworleyFB4  \
   nsdata ( "cworley"   , 0          , 0.647968 , 0.156623  , 0 ) // +/- 0.00587861
#define ns_falligator4 \
   nsdata ( "alligator" , 0          , 0.994222 , 0.117762  , 0 ) // +/- 0.0032
#define ns_valligator4 \
   nsdata ( "alligator" , 0          , 0.991346 , 0.125998  , 0 ) // +/- 0.0025
#define ns_fsine1      \
   nsdata ( "sine"      , 0          , 1        , 0.5       , 1 ) // +/- 0
#define ns_fsine2      \
   nsdata ( "sine"      , 0          , 1        , 0.5       , 1 ) // +/- 0
#define ns_fsine3      \
   nsdata ( "sine"      , 0          , 1        , 0.5       , 1 ) // +/- 0
#define ns_fsine4      \
   nsdata ( "sine"      , 0          , 1        , 0.5       , 1 ) // +/- 0

#define ns_vsine1      \
   nsdata ( "sine"      , 0          , 1        , 0.5       , 1 ) // +/- 0
#define ns_vsine2      \
   nsdata ( "sine"      , 0          , 1        , 0.5       , 1 ) // +/- 0
#define ns_vsine3      \
   nsdata ( "sine"      , 0          , 1        , 0.5       , 1 ) // +/- 0
#define ns_vsine4      \
   nsdata ( "sine"      , 0          , 1        , 0.5       , 1 ) // +/- 0

//------------------------------------------------------------------------------
// Pivots
// Not all noises are symetrical, so these functions return the natural 
// floor or "pivot" for each noise type.
//------------------------------------------------------------------------------
float npivot_symmetric(nsdata stats) { return stats.avg; }
float npivot_positive(nsdata stats)  { return 0.0125; }
float npivot_negative(nsdata stats)  { return 0.9875; }
#define npivot_sine      npivot_symmetric
#define npivot_perlin    npivot_symmetric
#define npivot_pperlin   npivot_symmetric
#define npivot_simplex   npivot_symmetric
#define npivot_sparse    npivot_symmetric
#define npivot_flow      npivot_symmetric
#define npivot_pflow     npivot_symmetric
#define npivot_worleyFA  npivot_positive
#define npivot_worleyFB  npivot_positive
#define npivot_mworleyFA npivot_positive
#define npivot_mworleyFB npivot_positive
#define npivot_cworleyFA npivot_positive
#define npivot_cworleyFB npivot_positive
#define npivot_alligator npivot_positive




// -------------------------------------------------------------------------------
// Wrappers
// These wrap the native noise functions, providing a single interface
// to all the noises used in pyro: They all take a vector4 position and period.
// They also stretch and transform the native output values to fit the [0,1]
// output range as tightly as possible.
// Naming convention: 
//       <f,v>nwrap_<p><name><dim>(vector4 position, period)
// Where:
//    -  The prefix 'f' denotes a function that returns a float, and 'v'
//       a vector.
//    -  A 'p' before the <name> denotes a periodic variant
//    -  The root <name> can be one of:
//       'perlin'    => noise()
//       'pperlin'   => pnoise()
//       'simplex'   => xnoise()
//       'sparse'    => snoise()
//       'flow'      => flownoise()
//       'pflow'     => flowpnoise()
//       'worleyFA'  => wnoise(), with function: f1
//       'worleyFB'  => wnoise(), with function: f2 - f1
//       'mworleyFA' => mwnoise(), with function: f1
//       'mworleyFB' => mwnoise(), with function: f2 - f1
//       'cworleyFA' => mwnoise(), with function: f1
//       'cworleyFB' => mwnoise(), with function: f2 - f1
//       'alligator' => anoise()
//    -  The dimension <dim> can be one of: 1, 2, 3, or 4
// -------------------------------------------------------------------------------

// Sine
//------------------------------------------------------------------------------
float  nwrap_sine(float p,per; float flow) { 
   return sin(p*C_PI)*.5+.5;
}
float  nwrap_sine(vector2 p,per; float flow) { 
   vector sv = sin(set(p.x,p.y,0)*C_PI)*.5+.5;
   return sv.x*sv.y;
}
float  nwrap_sine(vector p,per; float flow) { 
   return prod(sin(p*C_PI)*.5+.5);
}
float  nwrap_sine(vector4 p,per; float flow) { 
   return prod(sin(p*C_PI)*.5+.5);
}
vector nwrap_sine(float p,per; float flow) { 
   return set(nwrap_sine(p,per,flow),
              nwrap_sine(p-123.456,per,flow),
              nwrap_sine(p+123.456,per,flow));
}
vector nwrap_sine(vector2 p,per; float flow) { 
   return set(nwrap_sine(p,per,flow),
              nwrap_sine(p-123.456,per,flow),
              nwrap_sine(p+123.456,per,flow));
}
vector nwrap_sine(vector p,per; float flow) { 
    // for backwards compatibility
   vector4 pp = (vector4) p;
   pp.w=0;
   return set(nwrap_sine(pp,per,flow),
              nwrap_sine(pp-123.456,per,flow),
              nwrap_sine(pp+123.456,per,flow));
}
vector nwrap_sine(vector4 p,per; float flow) { 
   return set(nwrap_sine(p,per,flow),
              nwrap_sine(p-123.456,per,flow),
              nwrap_sine(p+123.456,per,flow));
}

// Perlin
//------------------------------------------------------------------------------
float  nwrap_perlin(float p,per; float flow) { 
   return fit(float(noise(p)), 
              ns_fperlin1.min, ns_fperlin1.max, 0,1);
}
float  nwrap_perlin(vector2 p,per; float flow) { 
   return fit(float(noise(p.x,p.y)), 
              ns_fperlin2.min, ns_fperlin2.max, 0,1);
}
float  nwrap_perlin(vector p,per; float flow) { 
   return fit(float(noise((vector)p)), 
              ns_fperlin3.min, ns_fperlin3.max, 0,1);
}
float  nwrap_perlin(vector4 p,per; float flow) { 
   return fit(float(noise(p)), 
              ns_fperlin4.min, ns_fperlin4.max, 0,1);
}
vector nwrap_perlin(float p,per; float flow) { 
   return fit(vector(noise(p)), 
              (vector)ns_vperlin1.min, 
              (vector)ns_vperlin1.max, VZERO,VONE);
}
vector nwrap_perlin(vector2 p,per; float flow) { 
   return fit(vector(noise(p.x,p.y)), 
              (vector)ns_vperlin2.min, 
              (vector)ns_vperlin2.max, VZERO,VONE);
}
vector nwrap_perlin(vector p,per; float flow) { 
   return fit(vector(noise((vector)p)), 
              (vector)ns_vperlin3.min, 
              (vector)ns_vperlin3.max, VZERO,VONE);
}
vector nwrap_perlin(vector4 p,per; float flow) { 
   return fit(vector(noise(p)), 
              (vector)ns_vperlin4.min, 
              (vector)ns_vperlin4.max, VZERO,VONE);
}


// Periodic Perlin
//------------------------------------------------------------------------------
float nwrap_pperlin(float p,per; float flow) { 
   return fit(float(pnoise(p,(int)per)), 
              ns_fpperlin1.min, ns_fpperlin1.max, 0,1);
}
float nwrap_pperlin(vector2 p,per; float flow) { 
   return fit(float(pnoise(p.x,p.y,(int)per.x,(int)per.y)),
              ns_fpperlin2.min, ns_fpperlin2.max, 0,1); 
}
float nwrap_pperlin(vector p,per; float flow) { 
   return fit(float(pnoise((vector)p,(vector)per)),
              ns_fpperlin3.min, ns_fpperlin3.max, 0,1);
}
float nwrap_pperlin(vector4 p,per; float flow) { 
   return fit(float(pnoise(p,per)),
              ns_fpperlin4.min, ns_fpperlin4.max, 0,1);
}
vector nwrap_pperlin(float p,per; float flow) { 
   return fit(vector(pnoise(p,(int)per)),
              (vector)ns_vpperlin1.min, 
              (vector)ns_vpperlin1.max, VZERO,VONE);
}
vector nwrap_pperlin(vector2 p,per; float flow) { 
   return fit(vector(pnoise(p.x,p.y,(int)per.x,(int)per.y)),
              (vector)ns_vpperlin2.min, 
              (vector)ns_vpperlin2.max, VZERO,VONE);
}
vector nwrap_pperlin(vector p,per; float flow) { 
   return fit(vector(pnoise((vector)p,(vector)per)),
              (vector)ns_vpperlin3.min, 
              (vector)ns_vpperlin3.max, VZERO,VONE);
}
vector nwrap_pperlin(vector4 p,per; float flow) { 
   return fit(vector(pnoise(p,per)),
              (vector)ns_vpperlin4.min, 
              (vector)ns_vpperlin4.max, VZERO,VONE);
}


// Simplex
//------------------------------------------------------------------------------
float  nwrap_simplex(float p,per; float flow) { 
   return fit(float(xnoise(p)),
              ns_fsimplex1.min, ns_fsimplex1.max, 0,1);
}
float  nwrap_simplex(vector2 p,per; float flow) { 
   return fit(float(xnoise(p.x,p.y)),
              ns_fsimplex2.min, ns_fsimplex2.max, 0,1);
}
float  nwrap_simplex(vector p,per; float flow) { 
   return fit(float(xnoise((vector)p)),
              ns_fsimplex3.min, ns_fsimplex3.max, 0,1);
}
float  nwrap_simplex(vector4 p,per; float flow) { 
   return fit(float(xnoise(p)),
              ns_fsimplex4.min, ns_fsimplex4.max, 0,1);
}
vector nwrap_simplex(float p,per; float flow) { 
   return fit(vector(xnoise(p)),
             (vector)ns_vsimplex1.min, 
             (vector)ns_vsimplex1.max, VZERO,VONE);
}
vector nwrap_simplex(vector2 p,per; float flow) { 
   return fit(vector(xnoise(p.x,p.y)),
             (vector)ns_vsimplex2.min, 
             (vector)ns_vsimplex2.max, VZERO,VONE);
}
vector nwrap_simplex(vector p,per; float flow) { 
   return fit(vector(xnoise(p)),
             (vector)ns_vsimplex3.min, 
             (vector)ns_vsimplex3.max, VZERO,VONE);
}
vector nwrap_simplex(vector4 p,per; float flow) { 
   return fit(vector(xnoise(p)),
             (vector)ns_vsimplex4.min, 
             (vector)ns_vsimplex4.max, VZERO,VONE);
}


// Sparse Convolution
//------------------------------------------------------------------------------
// Freq adjustment: scaled by 1.25 to match other noises better
float nwrap_sparse(float p,per; float flow) { 
   float chaos = float(snoise((vector)p*1.25));
   return fit(chaos,ns_fsparse1.min, ns_fsparse1.max, 0,1);
}
float nwrap_sparse(vector2 p,per; float flow) { 
   float chaos = float(snoise((vector)p*1.25));
   return fit(chaos,ns_fsparse2.min, ns_fsparse2.max, 0,1);
}
float nwrap_sparse(vector p,per; float flow) { 
   float chaos = float(snoise((vector)p*1.25));
   return fit(chaos,ns_fsparse3.min, ns_fsparse3.max, 0,1);
}
float nwrap_sparse(vector4 p,per; float flow) { 
   float chaos = float(snoise((vector)p*1.25));
   return fit(chaos,ns_fsparse4.min, ns_fsparse4.max, 0,1);
}
vector nwrap_sparse(float p,per; float flow) { 
   vector chaos = vector(snoise((vector)p*1.25));
   return fit(chaos,(vector)ns_vsparse1.min,
                    (vector)ns_vsparse1.max, VZERO,VONE);
}
vector nwrap_sparse(vector2 p,per; float flow) { 
   vector chaos = vector(snoise((vector)p*1.25));
   return fit(chaos,(vector)ns_vsparse2.min,
                    (vector)ns_vsparse2.max, VZERO,VONE);
}
vector nwrap_sparse(vector p,per; float flow) { 
   vector chaos = vector(snoise((vector)p*1.25));
   return fit(chaos,(vector)ns_vsparse3.min,
                    (vector)ns_vsparse3.max, VZERO,VONE);
}
vector nwrap_sparse(vector4 p,per; float flow) { 
   vector chaos = vector(snoise((vector)p*1.25));
   return fit(chaos,(vector)ns_vsparse4.min,
                    (vector)ns_vsparse4.max, VZERO,VONE);
}


// Flow
//------------------------------------------------------------------------------
float  nwrap_flow(float p,per; float flow) { 
   return fit(float(flownoise(p,0,flow)), 
              ns_fflow1.min, ns_fflow1.max, 0,1);
}
float  nwrap_flow(vector2 p,per; float flow) { 
   return fit(float(flownoise(p.x,p.y,flow)), 
              ns_fflow2.min, ns_fflow2.max, 0,1);
}
float  nwrap_flow(vector p,per; float flow) { 
   return fit(float(flownoise((vector)p,flow)), 
              ns_fflow3.min, ns_fflow3.max, 0,1);
}
float  nwrap_flow(vector4 p,per; float flow) { 
   return fit(float(flownoise(p,flow)), 
              ns_fflow4.min, ns_fflow4.max, 0,1);
}
vector nwrap_flow(float p,per; float flow) { 
   return fit(vector(flownoise(p,0,flow)),
             (vector)ns_vflow1.min,
             (vector)ns_vflow1.max, VZERO,VONE);
}
vector nwrap_flow(vector2 p,per; float flow) { 
   return fit(vector(flownoise(p.x,p.y,flow)),
             (vector)ns_vflow2.min,
             (vector)ns_vflow2.max, VZERO,VONE);
}
vector nwrap_flow(vector p,per; float flow) { 
   return fit(vector(flownoise((vector)p,flow)),
             (vector)ns_vflow3.min,
             (vector)ns_vflow3.max, VZERO,VONE);
}
vector nwrap_flow(vector4 p,per; float flow) { 
   return fit(vector(flownoise(p,flow)),
             (vector)ns_vflow4.min,
             (vector)ns_vflow4.max, VZERO,VONE);
}


// Periodic Flow
//------------------------------------------------------------------------------
float nwrap_pflow(float p,per; float flow) { 
   return fit(float(flowpnoise((vector)p,0,(int)per,0,flow)), 
              ns_fpflow1.min, ns_fpflow1.max, 0,1);
}
float nwrap_pflow(vector2 p,per; float flow) { 
   return fit(float(flowpnoise(p.x,p.y,(int)per.x,(int)per.y,flow)),
               ns_fpflow2.min, ns_fpflow2.max, 0,1); 
}
float nwrap_pflow(vector p,per; float flow) { 
   return fit(float(flowpnoise((vector)p,(vector)per,flow)),
               ns_fpflow3.min, ns_fpflow3.max, 0,1);
}
float nwrap_pflow(vector4 p,per; float flow) { 
   return fit(float(flowpnoise(p,per,flow)),
               ns_fpflow4.min, ns_fpflow4.max, 0,1);
}
vector nwrap_pflow(float p,per; float flow) { 
   return fit(vector(flowpnoise((vector)p,0,(int)per,0,flow)),
               (vector)ns_vpflow1.min,
               (vector)ns_vpflow1.max, VZERO,VONE);
}
vector nwrap_pflow(vector2 p,per; float flow) { 
   return fit(vector(flowpnoise(p.x,p.y,(int)per.x,(int)per.y,flow)),
               (vector)ns_vpflow2.min,
               (vector)ns_vpflow2.max, VZERO,VONE);
}
vector nwrap_pflow(vector p,per; float flow) { 
   return fit(vector(flowpnoise((vector)p,(vector)per,flow)),
               (vector)ns_vpflow3.min,
               (vector)ns_vpflow3.max, VZERO,VONE);
}
vector nwrap_pflow(vector4 p,per; float flow) { 
   return fit(vector(flowpnoise(p,per,flow)),
               (vector)ns_vpflow4.min,
               (vector)ns_vpflow4.max, VZERO,VONE);
}


// WorleyFA (f1)
//------------------------------------------------------------------------------
float nwrap_worleyFA(float p,per; float flow) { 
   int seed; float f1,f2; 
   wnoise(p,seed,f1,f2); 
   return fit(f1,ns_fworleyFA1.min,ns_fworleyFA1.max, 0,1);
}
float nwrap_worleyFA(vector2 p,per; float flow) { 
   int seed; float f1,f2; 
   wnoise(p.x,p.y,seed,f1,f2); 
   return fit(f1,ns_fworleyFA2.min,ns_fworleyFA2.max, 0,1);
}
float nwrap_worleyFA(vector p,per; float flow) { 
   int seed; float f1,f2; 
   wnoise((vector)p,seed,f1,f2);
   return fit(f1,ns_fworleyFA3.min,ns_fworleyFA3.max, 0,1);
}
float nwrap_worleyFA(vector4 p,per; float flow) { 
   int seed; float f1,f2; 
   wnoise(p,seed,f1,f2); 
   return fit(f1,ns_fworleyFA4.min,ns_fworleyFA4.max, 0,1);
}
vector nwrap_worleyFA(float p,per; float flow) {
   return set(nwrap_worleyFA(p,per,flow),
              nwrap_worleyFA(p+123.456,per,flow),
              nwrap_worleyFA(p-210.123,per,flow)); 
}
vector nwrap_worleyFA(vector2 p,per; float flow) {
   return set(nwrap_worleyFA(p,per,flow),
              nwrap_worleyFA(p+123.456,per,flow),
              nwrap_worleyFA(p-210.123,per,flow)); 
}
vector nwrap_worleyFA(vector p,per; float flow) {
   return set(nwrap_worleyFA(p,per,flow),
              nwrap_worleyFA(p+123.456,per,flow),
              nwrap_worleyFA(p-210.123,per,flow)); 
}
vector nwrap_worleyFA(vector4 p,per; float flow) {
   return set(nwrap_worleyFA(p,per,flow),
              nwrap_worleyFA(p+123.456,per,flow),
              nwrap_worleyFA(p-210.123,per,flow)); 
}


// WorleyFB (f2-f1)
//------------------------------------------------------------------------------
float nwrap_worleyFB(float p,per; float flow) { 
   int seed; float f1,f2; 
   wnoise(p,seed,f1,f2); 
   return fit(f2-f1, ns_fworleyFB1.min, ns_fworleyFB1.max, 0,1);
}
float nwrap_worleyFB(vector2 p,per; float flow) { 
   int seed; float f1,f2; 
   wnoise(p.x,p.y,seed,f1,f2); 
   return fit(f2-f1, ns_fworleyFB2.min, ns_fworleyFB2.max, 0,1);
}
float nwrap_worleyFB(vector p,per; float flow) { 
   int seed; float f1,f2; 
   wnoise((vector)p,seed,f1,f2);
   return fit(f2-f1, ns_fworleyFB3.min, ns_fworleyFB3.max, 0,1);
}
float nwrap_worleyFB(vector4 p,per; float flow) { 
   int seed; float f1,f2; 
   wnoise(p,seed,f1,f2);
   return fit(f2-f1, ns_fworleyFB4.min, ns_fworleyFB4.max, 0,1);
}
vector nwrap_worleyFB(float p,per; float flow) {
   return set(nwrap_worleyFB(p,per,flow),
              nwrap_worleyFB(p+123.456,per,flow),
              nwrap_worleyFB(p-210.123,per,flow)); 
}
vector nwrap_worleyFB(vector2 p,per; float flow) {
   return set(nwrap_worleyFB(p,per,flow),
              nwrap_worleyFB(p+123.456,per,flow),
              nwrap_worleyFB(p-210.123,per,flow)); 
}
vector nwrap_worleyFB(vector p,per; float flow) {
   return set(nwrap_worleyFB(p,per,flow),
              nwrap_worleyFB(p+123.456,per,flow),
              nwrap_worleyFB(p-210.123,per,flow)); 
}
vector nwrap_worleyFB(vector4 p,per; float flow) {
   return set(nwrap_worleyFB(p,per,flow),
              nwrap_worleyFB(p+123.456,per,flow),
              nwrap_worleyFB(p-210.123,per,flow)); 
}

// Manhattan WorleyFA (f1)
//------------------------------------------------------------------------------
float nwrap_mworleyFA(float p,per; float flow) { 
   int seed; float f1,f2; 
   mwnoise(p,seed,f1,f2); 
   return fit(f1,ns_fmworleyFA1.min,ns_fmworleyFA1.max, 0,1);
}
float nwrap_mworleyFA(vector2 p,per; float flow) { 
   int seed; float f1,f2; 
   mwnoise(p.x,p.y,seed,f1,f2); 
   return fit(f1,ns_fmworleyFA2.min,ns_fmworleyFA2.max, 0,1);
}
float nwrap_mworleyFA(vector p,per; float flow) { 
   int seed; float f1,f2; 
   mwnoise((vector)p,seed,f1,f2);
   return fit(f1,ns_fmworleyFA3.min,ns_fmworleyFA3.max, 0,1);
}
float nwrap_mworleyFA(vector4 p,per; float flow) { 
   int seed; float f1,f2; 
   mwnoise(p,seed,f1,f2); 
   return fit(f1,ns_fmworleyFA4.min,ns_fmworleyFA4.max, 0,1);
}
vector nwrap_mworleyFA(float p,per; float flow) {
   return set(nwrap_mworleyFA(p,per,flow),
              nwrap_mworleyFA(p+123.456,per,flow),
              nwrap_mworleyFA(p-210.123,per,flow)); 
}
vector nwrap_mworleyFA(vector2 p,per; float flow) {
   return set(nwrap_mworleyFA(p,per,flow),
              nwrap_mworleyFA(p+123.456,per,flow),
              nwrap_mworleyFA(p-210.123,per,flow)); 
}
vector nwrap_mworleyFA(vector p,per; float flow) {
   return set(nwrap_mworleyFA(p,per,flow),
              nwrap_mworleyFA(p+123.456,per,flow),
              nwrap_mworleyFA(p-210.123,per,flow)); 
}
vector nwrap_mworleyFA(vector4 p,per; float flow) {
   return set(nwrap_mworleyFA(p,per,flow),
              nwrap_mworleyFA(p+123.456,per,flow),
              nwrap_mworleyFA(p-210.123,per,flow)); 
}


// Manhattan WorleyFB (f2-f1)
//------------------------------------------------------------------------------
float nwrap_mworleyFB(float p,per; float flow) { 
   int seed; float f1,f2; 
   mwnoise(p,seed,f1,f2); 
   return fit(f2-f1, ns_fmworleyFB1.min, ns_fmworleyFB1.max, 0,1);
}
float nwrap_mworleyFB(vector2 p,per; float flow) { 
   int seed; float f1,f2; 
   mwnoise(p.x,p.y,seed,f1,f2); 
   return fit(f2-f1, ns_fmworleyFB2.min, ns_fmworleyFB2.max, 0,1);
}
float nwrap_mworleyFB(vector p,per; float flow) { 
   int seed; float f1,f2; 
   mwnoise((vector)p,seed,f1,f2);
   return fit(f2-f1, ns_fmworleyFB3.min, ns_fmworleyFB3.max, 0,1);
}
float nwrap_mworleyFB(vector4 p,per; float flow) { 
   int seed; float f1,f2; 
   mwnoise(p,seed,f1,f2);
   return fit(f2-f1, ns_fmworleyFB4.min, ns_fmworleyFB4.max, 0,1);
}
vector nwrap_mworleyFB(float p,per; float flow) {
   return set(nwrap_mworleyFB(p,per,flow),
              nwrap_mworleyFB(p+123.456,per,flow),
              nwrap_mworleyFB(p-210.123,per,flow)); 
}
vector nwrap_mworleyFB(vector2 p,per; float flow) {
   return set(nwrap_mworleyFB(p,per,flow),
              nwrap_mworleyFB(p+123.456,per,flow),
              nwrap_mworleyFB(p-210.123,per,flow)); 
}
vector nwrap_mworleyFB(vector p,per; float flow) {
   return set(nwrap_mworleyFB(p,per,flow),
              nwrap_mworleyFB(p+123.456,per,flow),
              nwrap_mworleyFB(p-210.123,per,flow)); 
}
vector nwrap_mworleyFB(vector4 p,per; float flow) {
   return set(nwrap_mworleyFB(p,per,flow),
              nwrap_mworleyFB(p+123.456,per,flow),
              nwrap_mworleyFB(p-210.123,per,flow)); 
}

// Chebyshev WorleyFA (f1)
//------------------------------------------------------------------------------
float nwrap_cworleyFA(float p,per; float flow) { 
   int seed; float f1,f2; 
   cwnoise(p,seed,f1,f2); 
   return fit(f1,ns_fcworleyFA1.min,ns_fcworleyFA1.max, 0,1);
}
float nwrap_cworleyFA(vector2 p,per; float flow) { 
   int seed; float f1,f2; 
   cwnoise(p.x,p.y,seed,f1,f2); 
   return fit(f1,ns_fcworleyFA2.min,ns_fcworleyFA2.max, 0,1);
}
float nwrap_cworleyFA(vector p,per; float flow) { 
   int seed; float f1,f2; 
   cwnoise((vector)p,seed,f1,f2);
   return fit(f1,ns_fcworleyFA3.min,ns_fcworleyFA3.max, 0,1);
}
float nwrap_cworleyFA(vector4 p,per; float flow) { 
   int seed; float f1,f2; 
   cwnoise(p,seed,f1,f2); 
   return fit(f1,ns_fcworleyFA4.min,ns_fcworleyFA4.max, 0,1);
}
vector nwrap_cworleyFA(float p,per; float flow) {
   return set(nwrap_cworleyFA(p,per,flow),
              nwrap_cworleyFA(p+123.456,per,flow),
              nwrap_cworleyFA(p-210.123,per,flow)); 
}
vector nwrap_cworleyFA(vector2 p,per; float flow) {
   return set(nwrap_cworleyFA(p,per,flow),
              nwrap_cworleyFA(p+123.456,per,flow),
              nwrap_cworleyFA(p-210.123,per,flow)); 
}
vector nwrap_cworleyFA(vector p,per; float flow) {
   return set(nwrap_cworleyFA(p,per,flow),
              nwrap_cworleyFA(p+123.456,per,flow),
              nwrap_cworleyFA(p-210.123,per,flow)); 
}
vector nwrap_cworleyFA(vector4 p,per; float flow) {
   return set(nwrap_cworleyFA(p,per,flow),
              nwrap_cworleyFA(p+123.456,per,flow),
              nwrap_cworleyFA(p-210.123,per,flow)); 
}

// Chebyshev WorleyFB (f2-f1)
//------------------------------------------------------------------------------
float nwrap_cworleyFB(float p,per; float flow) { 
   int seed; float f1,f2; 
   cwnoise(p,seed,f1,f2); 
   return fit(f2-f1, ns_fcworleyFB1.min, ns_fcworleyFB1.max, 0,1);
}
float nwrap_cworleyFB(vector2 p,per; float flow) { 
   int seed; float f1,f2; 
   cwnoise(p.x,p.y,seed,f1,f2); 
   return fit(f2-f1, ns_fcworleyFB2.min, ns_fcworleyFB2.max, 0,1);
}
float nwrap_cworleyFB(vector p,per; float flow) { 
   int seed; float f1,f2; 
   cwnoise((vector)p,seed,f1,f2);
   return fit(f2-f1, ns_fcworleyFB3.min, ns_fcworleyFB3.max, 0,1);
}
float nwrap_cworleyFB(vector4 p,per; float flow) { 
   int seed; float f1,f2; 
   cwnoise(p,seed,f1,f2);
   return fit(f2-f1, ns_fcworleyFB4.min, ns_fcworleyFB4.max, 0,1);
}
vector nwrap_cworleyFB(float p,per; float flow) {
   return set(nwrap_cworleyFB(p,per,flow),
              nwrap_cworleyFB(p+123.456,per,flow),
              nwrap_cworleyFB(p-210.123,per,flow)); 
}
vector nwrap_cworleyFB(vector2 p,per; float flow) {
   return set(nwrap_cworleyFB(p,per,flow),
              nwrap_cworleyFB(p+123.456,per,flow),
              nwrap_cworleyFB(p-210.123,per,flow)); 
}
vector nwrap_cworleyFB(vector p,per; float flow) {
   return set(nwrap_cworleyFB(p,per,flow),
              nwrap_cworleyFB(p+123.456,per,flow),
              nwrap_cworleyFB(p-210.123,per,flow)); 
}
vector nwrap_cworleyFB(vector4 p,per; float flow) {
   return set(nwrap_cworleyFB(p,per,flow),
              nwrap_cworleyFB(p+123.456,per,flow),
              nwrap_cworleyFB(p-210.123,per,flow)); 
}

// Alligator
//------------------------------------------------------------------------------
// Freq adjustment: scaled by 1.64 to better match the other noises
float nwrap_alligator(float p,per; float flow) { 
   return fit(float(anoise((vector)p*1.64)),
              ns_falligator1.min, ns_falligator1.max, 0,1);
}
float nwrap_alligator(vector2 p,per; float flow) { 
   return fit(float(anoise((vector)p*1.64)),
              ns_falligator2.min, ns_falligator2.max, 0,1);
}
float nwrap_alligator(vector p,per; float flow) { 
   return fit(float(anoise((vector)p*1.64)),
              ns_falligator3.min, ns_falligator3.max, 0,1);
}
float nwrap_alligator(vector4 p,per; float flow) { 
   return fit(float(anoise((vector)p*1.64)),
              ns_falligator4.min, ns_falligator4.max, 0,1);
}
vector nwrap_alligator(float p,per; float flow) { 
   return fit(vector(anoise((vector)p*1.64)),
              (vector)ns_valligator1.min,
              (vector)ns_valligator1.max, 0,1);
}
vector nwrap_alligator(vector2 p,per; float flow) { 
   return fit(vector(anoise((vector)p*1.64)),
              (vector)ns_valligator2.min,
              (vector)ns_valligator2.max, 0,1);
}
vector nwrap_alligator(vector p,per; float flow) { 
   return fit(vector(anoise((vector)p*1.64)),
              (vector)ns_valligator3.min,
              (vector)ns_valligator3.max, 0,1);
}
vector nwrap_alligator(vector4 p,per; float flow) { 
   return fit(vector(anoise((vector)p*1.64)),
              (vector)ns_valligator4.min,
              (vector)ns_valligator4.max, 0,1);
}



// Macro to construct a wrapper name
#define WNAME(name) CAT2(nwrap_,name)



//------------------------------------------------------------------------------
// Gradient Calculation -- by central differencing, which is expensive :/
// The only exception, at time of writing, is the new "simplex" noise, which
// has a gradient variant (and is therefore a lot faster than sampling the noise
// twice) -- RFE to add gradient calculation to all noise types!
//-----------------------------------------------------------------------------
#define GRADNAMEfloat(name)       CAT2(fngrad_,name)
#define GRADNAMEvector(name)      CAT2(vngrad_,name)
// dNoise/dP
#define NGRAD(rtype,ptype,name,dim) \
ptype GRADNAME##rtype(name) (ptype p,per; float flow,delta) { \
   float   d   = max(1e-5,abs(delta))*0.5; \
   float   di  = 1.0 / d; \
   int     i; \
   vector4 out=0,pd0=0,pd1=0; \
   for(i=0;i<dim;i++) { \
      pd0 = pd1 = p; pd0[i] += d; pd1[i] -= d; \
      rtype n0 = WNAME(name)(pd0,per,flow); \
      rtype n1 = WNAME(name)(pd1,per,flow); \
      out[i] = di * avg(n0 - n1); \
   } \
   return (ptype)out; \
}

// Sine
//-----------------------------------------------------------------------------
NGRAD(float ,float,  sine,1)
NGRAD(float ,vector2,sine,2)
NGRAD(float ,vector ,sine,3)
NGRAD(float ,vector4,sine,4)
NGRAD(vector,float,  sine,1)
NGRAD(vector,vector2,sine,2)
NGRAD(vector,vector ,sine,3)
NGRAD(vector,vector4,sine,4)


// Perlin
//-----------------------------------------------------------------------------
NGRAD(float ,float,  perlin,1)
NGRAD(float ,vector2,perlin,2)
NGRAD(float ,vector ,perlin,3)
NGRAD(float ,vector4,perlin,4)
NGRAD(vector,float,  perlin,1)
NGRAD(vector,vector2,perlin,2)
NGRAD(vector,vector ,perlin,3)
NGRAD(vector,vector4,perlin,4)


// Periodic Perlin
//-----------------------------------------------------------------------------
NGRAD(float ,float,  pperlin,1)
NGRAD(float ,vector2,pperlin,2)
NGRAD(float ,vector ,pperlin,3)
NGRAD(float ,vector4,pperlin,4)
NGRAD(vector,float,  pperlin,1)
NGRAD(vector,vector2,pperlin,2)
NGRAD(vector,vector ,pperlin,3)
NGRAD(vector,vector4,pperlin,4)


// Simplex
//-----------------------------------------------------------------------------
float fngrad_simplex(float p,per; float flow, delta) {
   float v1,g1;
   xnoised(p,v1,g1);
   return g1;
}
vector2 fngrad_simplex(vector2 p,per; float flow, delta) {
   float v1,g1,g2;
   xnoised(p.x,p.y,v1,g1,g2);
   return set(g1,g2);
}
vector fngrad_simplex(vector p,per; float flow, delta) {
   float v1,g1,g2,g3;
   xnoised((vector)p,v1,g1,g2,g3);
   return set(g1,g2,g3);
}
vector4 fngrad_simplex(vector4 p,per; float flow, delta) {
   float v1,g1,g2,g3,g4;
   xnoised(p,v1,g1,g2,g3,g4);
   return set(g1,g2,g3,g4);
}
float vngrad_simplex(float p,per; float flow, delta) {
   vector v1,g1;
   xnoised(p,v1,g1);
   return avg(g1);
}
vector2 vngrad_simplex(vector2 p,per; float flow, delta) {
   vector v1,g1,g2;
   xnoised(p.x,p.y,v1,g1,g2);
   return set(avg(g1),avg(g2));
}
vector vngrad_simplex(vector p,per; float flow, delta) {
   vector v1,g1,g2,g3;
   xnoised((vector)p,v1,g1,g2,g3);
   return set(avg(g1),avg(g2),avg(g3));
}
vector4 vngrad_simplex(vector4 p,per; float flow, delta) {
   vector v1,g1,g2,g3,g4;
   xnoised(p,v1,g1,g2,g3,g4);
   return set(avg(g1),avg(g2),avg(g3),avg(g4));
}

// Sparse Convolution
//-----------------------------------------------------------------------------
NGRAD(float ,float,  sparse,1)
NGRAD(float ,vector2,sparse,2)
NGRAD(float ,vector ,sparse,3)
NGRAD(float ,vector4,sparse,4)
NGRAD(vector,float,  sparse,1)
NGRAD(vector,vector2,sparse,2)
NGRAD(vector,vector ,sparse,3)
NGRAD(vector,vector4,sparse,4)


// Flow
//-----------------------------------------------------------------------------
NGRAD(float ,float,  flow,1)
NGRAD(float ,vector2,flow,2)
NGRAD(float ,vector ,flow,3)
NGRAD(float ,vector4,flow,4)
NGRAD(vector,float,  flow,1)
NGRAD(vector,vector2,flow,2)
NGRAD(vector,vector ,flow,3)
NGRAD(vector,vector4,flow,4)


// Flow Periodic
//-----------------------------------------------------------------------------
NGRAD(float ,float,  pflow,1)
NGRAD(float ,vector2,pflow,2)
NGRAD(float ,vector ,pflow,3)
NGRAD(float ,vector4,pflow,4)
NGRAD(vector,float,  pflow,1)
NGRAD(vector,vector2,pflow,2)
NGRAD(vector,vector ,pflow,3)
NGRAD(vector,vector4,pflow,4)



// Cellular (Worley) -- "FA" = f1
//-----------------------------------------------------------------------------
NGRAD(float ,float,  worleyFA,1)
NGRAD(float ,vector2,worleyFA,2)
NGRAD(float ,vector ,worleyFA,3)
NGRAD(float ,vector4,worleyFA,4)
NGRAD(vector,float,  worleyFA,1)
NGRAD(vector,vector2,worleyFA,2)
NGRAD(vector,vector ,worleyFA,3)
NGRAD(vector,vector4,worleyFA,4)


// Cellular (Worley) -- "FB" = f2-f1
//-----------------------------------------------------------------------------
NGRAD(float ,float,  worleyFB,1)
NGRAD(float ,vector2,worleyFB,2)
NGRAD(float ,vector ,worleyFB,3)
NGRAD(float ,vector4,worleyFB,4)
NGRAD(vector,float,  worleyFB,1)
NGRAD(vector,vector2,worleyFB,2)
NGRAD(vector,vector ,worleyFB,3)
NGRAD(vector,vector4,worleyFB,4)


// Cellular (Worley) Manhattan metric-- "FA" = f1
//-----------------------------------------------------------------------------
NGRAD(float ,float,  mworleyFA,1)
NGRAD(float ,vector2,mworleyFA,2)
NGRAD(float ,vector ,mworleyFA,3)
NGRAD(float ,vector4,mworleyFA,4)
NGRAD(vector,float,  mworleyFA,1)
NGRAD(vector,vector2,mworleyFA,2)
NGRAD(vector,vector ,mworleyFA,3)
NGRAD(vector,vector4,mworleyFA,4)


// Cellular (Worley) Manhattan metric-- "FB" = f2-f1
//-----------------------------------------------------------------------------
NGRAD(float ,float,  mworleyFB,1)
NGRAD(float ,vector2,mworleyFB,2)
NGRAD(float ,vector ,mworleyFB,3)
NGRAD(float ,vector4,mworleyFB,4)
NGRAD(vector,float,  mworleyFB,1)
NGRAD(vector,vector2,mworleyFB,2)
NGRAD(vector,vector ,mworleyFB,3)
NGRAD(vector,vector4,mworleyFB,4)


// Cellular (Worley) Chebyshev metric-- "FA" = f1
//-----------------------------------------------------------------------------
NGRAD(float ,float,  cworleyFA,1)
NGRAD(float ,vector2,cworleyFA,2)
NGRAD(float ,vector ,cworleyFA,3)
NGRAD(float ,vector4,cworleyFA,4)
NGRAD(vector,float,  cworleyFA,1)
NGRAD(vector,vector2,cworleyFA,2)
NGRAD(vector,vector ,cworleyFA,3)
NGRAD(vector,vector4,cworleyFA,4)


// Cellular (Worley) Chebyshev metric-- "FB" = f2-f1
//-----------------------------------------------------------------------------
NGRAD(float ,float,  cworleyFB,1)
NGRAD(float ,vector2,cworleyFB,2)
NGRAD(float ,vector ,cworleyFB,3)
NGRAD(float ,vector4,cworleyFB,4)
NGRAD(vector,float,  cworleyFB,1)
NGRAD(vector,vector2,cworleyFB,2)
NGRAD(vector,vector ,cworleyFB,3)
NGRAD(vector,vector4,cworleyFB,4)


// Alligator
//-----------------------------------------------------------------------------
NGRAD(float ,float,  alligator,1)
NGRAD(float ,vector2,alligator,2)
NGRAD(float ,vector ,alligator,3)
NGRAD(float ,vector4,alligator,4)
NGRAD(vector,float,  alligator,1)
NGRAD(vector,vector2,alligator,2)
NGRAD(vector,vector ,alligator,3)
NGRAD(vector,vector4,alligator,4)

#undef NGRAD



// Common exports for all bases
// x_avg: the calculated median for the current parameterization
// x_off: the sample's positional offset due to warps (in current space)
#define BXPARMS export float  x_avg; \
                export vector x_off
#define BXARGS  x_avg, x_off

// Common exports for all fractals
// x_oct: actual number of octaves computed for this sample
#define FXPARMS export float  x_oct
#define FXARGS  x_oct


// Common parms/args for all bases
#define BPARMS(FW, PTYPE) int inv,fold,accl,accg; float FW,expon; PTYPE per; \
                   int dolw,dogw; float disp,dfreq,gflow,flow
#define BARGS(FW)  inv,fold,accl,accg,FW,expon,per, \
                   dolw,dogw,disp,dfreq,gflow,flow

// Common parms/args for all fractals
#define FPARMS float maxoctaves,lacunarity,gain
#define FARGS	     maxoctaves,lacunarity,gain



#define FADE(rt,navg,f,fw) \
   ( rt(lerp(f,(rt)navg,smooth(0.35, 0.85, fw))) )

#define FADEF(rt,navg,f,fw) \
   ( rt(lerp(rt(f),(rt)navg,smooth(0.35, 0.85, fw))) )



//-----------------------------------------------------------------------------
// Basic noise correction (complement, fold, gamma)
//-----------------------------------------------------------------------------
float noise_correct(float val; nsdata stats; BPARMS(fw,vector4)) {
   if(expon==0) return 1;
   float navg = stats.avg;
   float out  = val;
   if(fold && stats.symmetric) out = abs(out-navg) / max(navg,1-navg);
   if(inv)      out = 1 - out;
   if(expon!=1) out = pow(out,expon);
   return out;
}
vector noise_correct(vector val; nsdata stats; BPARMS(fw,vector4)) {
   if(expon==0) return 1;
   float  navg = stats.avg;
   vector out  = val;
   if(fold && stats.symmetric) out = abs(out-navg) / max(navg,1-navg);
   if(inv)      out = VONE - out;
   if(expon!=1) out = pow(out,expon);
   return out;
}
float navg_correct(nsdata stats; BPARMS(fw,vector4)) {
   if(expon==0) return 1;
   float navg = stats.avg;
   if(fold && stats.symmetric) navg = navg*navg;
   if(inv)      navg = 1 - navg;
   if(expon!=1) navg = pow(navg,expon);
   return navg;
}
float navg_correct(float tgt; nsdata stats; BPARMS(fw,vector4)) {
   if(expon==0) return 1;
   float navg = tgt;
   if(fold && stats.symmetric) navg = navg*navg;
   if(inv)      navg = 1 - navg;
   if(expon!=1) navg = pow(navg,expon);
   return navg;
}




//-----------------------------------------------------------------------------
// Noise Bases
// These are the base functions called by the fractal generators and constitute
// the front-end API to the pyro noises. They extend the low level wrappers by 
// adding: 
// 1) basic range correction (complement,fold,exponent) 
// 2) lattice warp (noise distortion of the lookup space)
// 3) gradient warp (lookup offset along gradient of current sample)
//-----------------------------------------------------------------------------

#define BNAME(name)		   CAT2(noise_,name)
#define SNAMEfloat(name,dim)       CAT3(ns_f,name,dim)
#define SNAMEvector(name,dim)      CAT3(ns_v,name,dim)

// The basis definition
//-----------------------------------------------------------------------------
#define DEF_BASE(rtype,ptype,name,dim) \
rtype BNAME(name) (ptype pp; BPARMS(fw,ptype); BXPARMS ) { \
   if(expon==0) { x_avg=1; x_off=0; return 1;} \
   nsdata stats = SNAME##rtype(name,dim); \
   x_off = 0; \
   float med = x_avg = stats.avg; \
   ptype p = pp; \
   ptype dpl=0,dpg=0; \
\
   if(disp*dolw!=0) { \
      ptype pl = (pp+123.456)*dfreq; \
      rtype ln = WNAME(name)(pl,per,flow); \
      dpl = (ptype)(disp*(ln - med) / \
            max(med,1-med)); \
   }  \
   if(gflow*dogw!=0) { \
      dpg = gflow * GRADNAME##rtype(name) (pp,per,flow,0.05); \
   }  \
   pp += accl*dpl + accg*dpg; \
   x_off = (vector)(dpl + dpg); \
   x_avg = med; \
\
   p  += (ptype)x_off; \
   rtype out = WNAME(name)(p,per,flow); \
   out   = noise_correct(out,stats,BARGS(fw)); \
   x_avg = navg_correct(stats,BARGS(fw)); \
   float fw2 = fw*(abs(disp*dfreq)+1); \
   out = FADE(rtype,x_avg,out,fw2); \
   return out; \
}

// Sine
//-----------------------------------------------------------------------------
DEF_BASE(float ,float,  sine,1)
DEF_BASE(float ,vector2,sine,2)
DEF_BASE(float ,vector ,sine,3)
DEF_BASE(float ,vector4,sine,4)
DEF_BASE(vector,float,  sine,1)
DEF_BASE(vector,vector2,sine,2)
DEF_BASE(vector,vector ,sine,3)
DEF_BASE(vector,vector4,sine,4)


// Perlin
//-----------------------------------------------------------------------------
DEF_BASE(float ,float,  perlin,1)
DEF_BASE(float ,vector2,perlin,2)
DEF_BASE(float ,vector ,perlin,3)
DEF_BASE(float ,vector4,perlin,4)
DEF_BASE(vector,float,  perlin,1)
DEF_BASE(vector,vector2,perlin,2)
DEF_BASE(vector,vector ,perlin,3)
DEF_BASE(vector,vector4,perlin,4)


// Periodic Perlin
//-----------------------------------------------------------------------------
DEF_BASE(float ,float,  pperlin,1)
DEF_BASE(float ,vector2,pperlin,2)
DEF_BASE(float ,vector ,pperlin,3)
DEF_BASE(float ,vector4,pperlin,4)
DEF_BASE(vector,float,  pperlin,1)
DEF_BASE(vector,vector2,pperlin,2)
DEF_BASE(vector,vector ,pperlin,3)
DEF_BASE(vector,vector4,pperlin,4)


// Simplex
//-----------------------------------------------------------------------------
DEF_BASE(float ,float,  simplex,1)
DEF_BASE(float ,vector2,simplex,2)
DEF_BASE(float ,vector ,simplex,3)
DEF_BASE(float ,vector4,simplex,4)
DEF_BASE(vector,float,  simplex,1)
DEF_BASE(vector,vector2,simplex,2)
DEF_BASE(vector,vector ,simplex,3)
DEF_BASE(vector,vector4,simplex,4)


// Sparse Convolution
//-----------------------------------------------------------------------------
DEF_BASE(float ,float,  sparse,1)
DEF_BASE(float ,vector2,sparse,2)
DEF_BASE(float ,vector ,sparse,3)
DEF_BASE(float ,vector4,sparse,4)
DEF_BASE(vector,float,  sparse,1)
DEF_BASE(vector,vector2,sparse,2)
DEF_BASE(vector,vector ,sparse,3)
DEF_BASE(vector,vector4,sparse,4)


// Flow
//-----------------------------------------------------------------------------
DEF_BASE(float ,float,  flow,1)
DEF_BASE(float ,vector2,flow,2)
DEF_BASE(float ,vector ,flow,3)
DEF_BASE(float ,vector4,flow,4)
DEF_BASE(vector,float,  flow,1)
DEF_BASE(vector,vector2,flow,2)
DEF_BASE(vector,vector ,flow,3)
DEF_BASE(vector,vector4,flow,4)


// Flow Periodic
//-----------------------------------------------------------------------------
DEF_BASE(float ,float,  pflow,1)
DEF_BASE(float ,vector2,pflow,2)
DEF_BASE(float ,vector ,pflow,3)
DEF_BASE(float ,vector4,pflow,4)
DEF_BASE(vector,float,  pflow,1)
DEF_BASE(vector,vector2,pflow,2)
DEF_BASE(vector,vector ,pflow,3)
DEF_BASE(vector,vector4,pflow,4)



// Cellular (Worley) -- "FA" = f1
//-----------------------------------------------------------------------------
DEF_BASE(float ,float,  worleyFA,1)
DEF_BASE(float ,vector2,worleyFA,2)
DEF_BASE(float ,vector ,worleyFA,3)
DEF_BASE(float ,vector4,worleyFA,4)
DEF_BASE(vector,float,  worleyFA,1)
DEF_BASE(vector,vector2,worleyFA,2)
DEF_BASE(vector,vector ,worleyFA,3)
DEF_BASE(vector,vector4,worleyFA,4)


// Cellular (Worley) -- "FB" = f2-f1
//-----------------------------------------------------------------------------
DEF_BASE(float ,float,  worleyFB,1)
DEF_BASE(float ,vector2,worleyFB,2)
DEF_BASE(float ,vector ,worleyFB,3)
DEF_BASE(float ,vector4,worleyFB,4)
DEF_BASE(vector,float,  worleyFB,1)
DEF_BASE(vector,vector2,worleyFB,2)
DEF_BASE(vector,vector ,worleyFB,3)
DEF_BASE(vector,vector4,worleyFB,4)


// Cellular (Worley) Manhattan metric-- "FA" = f1
//-----------------------------------------------------------------------------
DEF_BASE(float ,float,  mworleyFA,1)
DEF_BASE(float ,vector2,mworleyFA,2)
DEF_BASE(float ,vector ,mworleyFA,3)
DEF_BASE(float ,vector4,mworleyFA,4)
DEF_BASE(vector,float,  mworleyFA,1)
DEF_BASE(vector,vector2,mworleyFA,2)
DEF_BASE(vector,vector ,mworleyFA,3)
DEF_BASE(vector,vector4,mworleyFA,4)


// Cellular (Worley) Manhattan metric-- "FB" = f2-f1
//-----------------------------------------------------------------------------
DEF_BASE(float ,float,  mworleyFB,1)
DEF_BASE(float ,vector2,mworleyFB,2)
DEF_BASE(float ,vector ,mworleyFB,3)
DEF_BASE(float ,vector4,mworleyFB,4)
DEF_BASE(vector,float,  mworleyFB,1)
DEF_BASE(vector,vector2,mworleyFB,2)
DEF_BASE(vector,vector ,mworleyFB,3)
DEF_BASE(vector,vector4,mworleyFB,4)


// Cellular (Worley) Chebyshev metric-- "FA" = f1
//-----------------------------------------------------------------------------
DEF_BASE(float ,float,  cworleyFA,1)
DEF_BASE(float ,vector2,cworleyFA,2)
DEF_BASE(float ,vector ,cworleyFA,3)
DEF_BASE(float ,vector4,cworleyFA,4)
DEF_BASE(vector,float,  cworleyFA,1)
DEF_BASE(vector,vector2,cworleyFA,2)
DEF_BASE(vector,vector ,cworleyFA,3)
DEF_BASE(vector,vector4,cworleyFA,4)


// Cellular (Worley) Chebyshev metric-- "FB" = f2-f1
//-----------------------------------------------------------------------------
DEF_BASE(float ,float,  cworleyFB,1)
DEF_BASE(float ,vector2,cworleyFB,2)
DEF_BASE(float ,vector ,cworleyFB,3)
DEF_BASE(float ,vector4,cworleyFB,4)
DEF_BASE(vector,float,  cworleyFB,1)
DEF_BASE(vector,vector2,cworleyFB,2)
DEF_BASE(vector,vector ,cworleyFB,3)
DEF_BASE(vector,vector4,cworleyFB,4)


// Alligator
//-----------------------------------------------------------------------------
DEF_BASE(float ,float,  alligator,1)
DEF_BASE(float ,vector2,alligator,2)
DEF_BASE(float ,vector ,alligator,3)
DEF_BASE(float ,vector4,alligator,4)
DEF_BASE(vector,float,  alligator,1)
DEF_BASE(vector,vector2,alligator,2)
DEF_BASE(vector,vector ,alligator,3)
DEF_BASE(vector,vector4,alligator,4)




//-----------------------------------------------------------------------------
// Normalization Parameters
// The main normalization is geared toward noises with gaussian distributions,
// like perlin, but needs to be adjusted for other types with significantly
// different distributions.
//-----------------------------------------------------------------------------
struct nndata {
   float scale;
   float blend;
}

#define nn_perlin    nndata(1.6, 1.8)
#define nn_pperlin   nn_perlin
#define nn_simplex   nn_perlin
#define nn_sparse    nn_perlin
#define nn_flow      nn_perlin
#define nn_pflow     nn_perlin
#define nn_worleyFA  nn_perlin
#define nn_worleyFB  nn_perlin
#define nn_mworleyFA nn_perlin
#define nn_mworleyFB nn_perlin
#define nn_cworleyFA nn_perlin
#define nn_cworleyFB nn_perlin
#define nn_alligator nn_perlin
#define nn_sine      nndata(6.0, 2.5)



//-----------------------------------------------------------------------------
// Attenuation at Extrema
// Very rarely -- mostly caused by extreme exponentiaion and/or lacunarity<=1,
// the normalization will be underestimated slightly.
// These functions are used to add a small amount of compression to the very
// extremes of the fractalized signal to bring in those pesky outliers.
//-----------------------------------------------------------------------------

float nattencomp(float x,k,e) {
   return pow(1 - x*x, e*.5) * (1-k) + k;
}
float natten(float n, k,e) {
   float x = n * 2 - 1;
   return nattencomp(x,k,e);
}
vector natten(vector n; float k,e) {
   vector x = n * 2 - 1;
   return set(nattencomp(x.x,k,e),
              nattencomp(x.y,k,e),
              nattencomp(x.z,k,e));
}
vector4 natten(vector n; float k,e) {
   vector4 x = n * (vector4)2 - (vector4)1;
   return set(nattencomp(x.x,k,e),
              nattencomp(x.y,k,e),
              nattencomp(x.z,k,e),
              nattencomp(x.w,k,e));
}








// -------------------------------------------------------------------------------
// Fractal Constructions
// Naming: <frac_type>_<base_name>
//     eg: fBm_perlin()
// -------------------------------------------------------------------------------
#define FNAME(fractal,name)  CAT2(fractal##_,name)


// fBm
// -----------------------------------------------
#define FRAC_FBM(RTYPE,PTYPE,NAME,DIM) \
RTYPE FNAME(fBm,NAME) (PTYPE p; BPARMS(fsize,PTYPE); BXPARMS; \
                  FPARMS; FXPARMS ) \
{ \
   nsdata stats = SNAME##RTYPE(NAME,DIM); \
   nndata norm  = nn_##NAME; \
   string name  = stats.name; \
   float  oct   = x_oct = 0; \
   float  med   = x_avg = navg_correct(stats,BARGS(fsize)); \
   float  pvt   = navg_correct(npivot_##NAME(stats),stats,BARGS(fsize)); \
   RTYPE  base  = BNAME(NAME) (p,BARGS(fsize),BXARGS); \
   if(fsize>=1 || oct>maxoctaves || lacunarity*maxoctaves*gain==0) \
      return base; \
\
   PTYPE   pp    = p; \
   RTYPE   sbase = base - pvt; \
   float   lac   = abs(lacunarity); \
   float   g     = gain * min(1,lac); \
   float   fw    = fsize*(1+fold), fwold;\
   RTYPE   out   = sbase, chaos; \
   RTYPE   wsum  = 1, w = 1; \
   RTYPE   wk    = lerp(10, fold*stats.symmetric ? 10 : norm.scale, \
                        smooth(1,norm.blend,lac)); \
\
   do { \
      w *= g; \
      oct += 1; \
      if(oct>maxoctaves) { \
         oct = maxoctaves; \
         w  *= bias(frac(maxoctaves),0.5-0.05); \
      } \
      fwold = fw; pp *= lac; fw *= lac; \
      if(fw>1) { \
         float k = (1-fwold) / (fw-fwold); \
         oct -= 1-k; \
         w   *= k; \
      } \
      chaos = BNAME(NAME) (pp,BARGS(fsize),BXARGS); \
      chaos -= pvt; \
      out  += w * chaos; \
      wsum += wk * w*w; \
\
   } while(oct<maxoctaves && fw<1.0); \
\
   x_oct = oct; \
   if(wsum!=0) out /= sqrt(wsum); \
   out = out*fit(g,0,0.1,1,natten(out+pvt,0.9,2)) + pvt; \
   out = clamp(out,(RTYPE)0,(RTYPE)1); /* safety - there should be no clipping */\
   return out; \
}


// Definitions for: fBm, 1D
//------------------------------------------------------------------------------
FRAC_FBM( float, float, sine,      1)
FRAC_FBM( float, float, perlin,    1)
FRAC_FBM( float, float, pperlin,   1)
FRAC_FBM( float, float, simplex,   1)
FRAC_FBM( float, float, sparse,    1)
FRAC_FBM( float, float, flow,      1)
FRAC_FBM( float, float, pflow,     1)
FRAC_FBM( float, float, worleyFA,  1)
FRAC_FBM( float, float, worleyFB,  1)
FRAC_FBM( float, float, mworleyFA, 1)
FRAC_FBM( float, float, mworleyFB, 1)
FRAC_FBM( float, float, cworleyFA, 1)
FRAC_FBM( float, float, cworleyFB, 1)
FRAC_FBM( float, float, alligator, 1)
FRAC_FBM(vector, float, sine,      1)
FRAC_FBM(vector, float, perlin,    1)
FRAC_FBM(vector, float, pperlin,   1)
FRAC_FBM(vector, float, simplex,   1)
FRAC_FBM(vector, float, sparse,    1)
FRAC_FBM(vector, float, flow,      1)
FRAC_FBM(vector, float, pflow,     1)
FRAC_FBM(vector, float, worleyFA,  1)
FRAC_FBM(vector, float, worleyFB,  1)
FRAC_FBM(vector, float, mworleyFA, 1)
FRAC_FBM(vector, float, mworleyFB, 1)
FRAC_FBM(vector, float, cworleyFA, 1)
FRAC_FBM(vector, float, cworleyFB, 1)
FRAC_FBM(vector, float, alligator, 1)

// Definitions for: fBm, 2D
//------------------------------------------------------------------------------
FRAC_FBM( float, vector2, sine,      2)
FRAC_FBM( float, vector2, perlin,    2)
FRAC_FBM( float, vector2, pperlin,   2)
FRAC_FBM( float, vector2, simplex,   2)
FRAC_FBM( float, vector2, sparse,    2)
FRAC_FBM( float, vector2, flow,      2)
FRAC_FBM( float, vector2, pflow,     2)
FRAC_FBM( float, vector2, worleyFA,  2)
FRAC_FBM( float, vector2, worleyFB,  2)
FRAC_FBM( float, vector2, mworleyFA, 2)
FRAC_FBM( float, vector2, mworleyFB, 2)
FRAC_FBM( float, vector2, cworleyFA, 2)
FRAC_FBM( float, vector2, cworleyFB, 2)
FRAC_FBM( float, vector2, alligator, 2)
FRAC_FBM(vector, vector2, sine,      2)
FRAC_FBM(vector, vector2, perlin,    2)
FRAC_FBM(vector, vector2, pperlin,   2)
FRAC_FBM(vector, vector2, simplex,   2)
FRAC_FBM(vector, vector2, sparse,    2)
FRAC_FBM(vector, vector2, flow,      2)
FRAC_FBM(vector, vector2, pflow,     2)
FRAC_FBM(vector, vector2, worleyFA,  2)
FRAC_FBM(vector, vector2, worleyFB,  2)
FRAC_FBM(vector, vector2, mworleyFA, 2)
FRAC_FBM(vector, vector2, mworleyFB, 2)
FRAC_FBM(vector, vector2, cworleyFA, 2)
FRAC_FBM(vector, vector2, cworleyFB, 2)
FRAC_FBM(vector, vector2, alligator, 2)

// Definitions for: fBm, 3D
//------------------------------------------------------------------------------
FRAC_FBM( float, vector, sine,      3)
FRAC_FBM( float, vector, perlin,    3)
FRAC_FBM( float, vector, pperlin,   3)
FRAC_FBM( float, vector, simplex,   3)
FRAC_FBM( float, vector, sparse,    3)
FRAC_FBM( float, vector, flow,      3)
FRAC_FBM( float, vector, pflow,     3)
FRAC_FBM( float, vector, worleyFA,  3)
FRAC_FBM( float, vector, worleyFB,  3)
FRAC_FBM( float, vector, mworleyFA, 3)
FRAC_FBM( float, vector, mworleyFB, 3)
FRAC_FBM( float, vector, cworleyFA, 3)
FRAC_FBM( float, vector, cworleyFB, 3)
FRAC_FBM( float, vector, alligator, 3)
FRAC_FBM(vector, vector, sine,      3)
FRAC_FBM(vector, vector, perlin,    3)
FRAC_FBM(vector, vector, pperlin,   3)
FRAC_FBM(vector, vector, simplex,   3)
FRAC_FBM(vector, vector, sparse,    3)
FRAC_FBM(vector, vector, flow,      3)
FRAC_FBM(vector, vector, pflow,     3)
FRAC_FBM(vector, vector, worleyFA,  3)
FRAC_FBM(vector, vector, worleyFB,  3)
FRAC_FBM(vector, vector, mworleyFA, 3)
FRAC_FBM(vector, vector, mworleyFB, 3)
FRAC_FBM(vector, vector, cworleyFA, 3)
FRAC_FBM(vector, vector, cworleyFB, 3)
FRAC_FBM(vector, vector, alligator, 3)

// Definitions for: fBm, 4D
//------------------------------------------------------------------------------
FRAC_FBM( float, vector4, sine,      4)
FRAC_FBM( float, vector4, perlin,    4)
FRAC_FBM( float, vector4, pperlin,   4)
FRAC_FBM( float, vector4, simplex,   4)
FRAC_FBM( float, vector4, sparse,    4)
FRAC_FBM( float, vector4, flow,      4)
FRAC_FBM( float, vector4, pflow,     4)
FRAC_FBM( float, vector4, worleyFA,  4)
FRAC_FBM( float, vector4, worleyFB,  4)
FRAC_FBM( float, vector4, mworleyFA, 4)
FRAC_FBM( float, vector4, mworleyFB, 4)
FRAC_FBM( float, vector4, cworleyFA, 4)
FRAC_FBM( float, vector4, cworleyFB, 4)
FRAC_FBM( float, vector4, alligator, 4)
FRAC_FBM(vector, vector4, sine,      4)
FRAC_FBM(vector, vector4, perlin,    4)
FRAC_FBM(vector, vector4, pperlin,   4)
FRAC_FBM(vector, vector4, simplex,   4)
FRAC_FBM(vector, vector4, sparse,    4)
FRAC_FBM(vector, vector4, flow,      4)
FRAC_FBM(vector, vector4, pflow,     4)
FRAC_FBM(vector, vector4, worleyFA,  4)
FRAC_FBM(vector, vector4, worleyFB,  4)
FRAC_FBM(vector, vector4, mworleyFA, 4)
FRAC_FBM(vector, vector4, mworleyFB, 4)
FRAC_FBM(vector, vector4, cworleyFA, 4)
FRAC_FBM(vector, vector4, cworleyFB, 4)
FRAC_FBM(vector, vector4, alligator, 4)
#undef FRAC_FBM




// -----------------------------------------------
// Multifractal Terrain (mfT)
// Successive frequencies are scaled by the 
// function's local value
// -----------------------------------------------
#define FRAC_MFT(RTYPE,PTYPE,NAME,DIM) \
RTYPE FNAME(mfT,NAME) (PTYPE p; BPARMS(fsize,PTYPE); BXPARMS; \
                  FPARMS; FXPARMS ) \
{ \
   nsdata stats = SNAME##RTYPE(NAME,DIM); \
   nndata norm  = nn_##NAME; \
   string name  = stats.name; \
   float  oct   = x_oct = 0; \
   float  med   = x_avg = navg_correct(stats,BARGS(fsize)); \
   float  pvt   = navg_correct(npivot_##NAME(stats),stats,BARGS(fsize)); \
   RTYPE  base  = BNAME(NAME) (p,BARGS(fsize),BXARGS); \
   if(fsize>=1 || oct>maxoctaves || lacunarity*maxoctaves*gain==0) \
      return base; \
\
   PTYPE   pp    = p; \
   RTYPE   sbase = base - pvt; \
   float   lac   = abs(lacunarity); \
   float   g     = abs(gain) * min(1,lac); \
   float   fw    = fsize*(1+fold), fwold;\
   RTYPE   out   = sbase, chaos; \
   RTYPE   wsum  = 1, w = base; \
   RTYPE   wk    = lerp(10, fold*stats.symmetric ? 10 : norm.scale, \
                        smooth(1,norm.blend,lac)); \
\
   do { \
      w *= g; \
      oct += 1; \
      if(oct>maxoctaves) { \
         oct = maxoctaves; \
         w  *= bias(frac(maxoctaves),0.5-0.05); \
      } \
      fwold = fw; pp *= lac; fw *= lac; \
      if(fw>1) { \
         float k = (1-fwold) / (fw-fwold); \
         oct -= 1-k; \
         w   *= k; \
      } \
      chaos = BNAME(NAME) (pp,BARGS(fsize),BXARGS); \
      chaos -= pvt; \
      out  += w * chaos; \
      wsum += wk * w*w; \
\
   } while(oct<maxoctaves && fw<1.0); \
\
   x_oct = oct; \
   if(wsum!=0) out /= sqrt(wsum); \
   out = out*fit(g,0,0.1,1,natten(out+pvt,0.9,2)) + pvt; \
   out = clamp(out,(RTYPE)0,(RTYPE)1); /* safety - there should be no clipping */\
   return out; \
}


// Definitions for: mfT, 1D
//------------------------------------------------------------------------------
FRAC_MFT( float, float, sine,      1)
FRAC_MFT( float, float, perlin,    1)
FRAC_MFT( float, float, pperlin,   1)
FRAC_MFT( float, float, simplex,   1)
FRAC_MFT( float, float, sparse,    1)
FRAC_MFT( float, float, flow,      1)
FRAC_MFT( float, float, pflow,     1)
FRAC_MFT( float, float, worleyFA,  1)
FRAC_MFT( float, float, worleyFB,  1)
FRAC_MFT( float, float, mworleyFA, 1)
FRAC_MFT( float, float, mworleyFB, 1)
FRAC_MFT( float, float, cworleyFA, 1)
FRAC_MFT( float, float, cworleyFB, 1)
FRAC_MFT( float, float, alligator, 1)
FRAC_MFT(vector, float, sine,      1)
FRAC_MFT(vector, float, perlin,    1)
FRAC_MFT(vector, float, pperlin,   1)
FRAC_MFT(vector, float, simplex,   1)
FRAC_MFT(vector, float, sparse,    1)
FRAC_MFT(vector, float, flow,      1)
FRAC_MFT(vector, float, pflow,     1)
FRAC_MFT(vector, float, worleyFA,  1)
FRAC_MFT(vector, float, worleyFB,  1)
FRAC_MFT(vector, float, mworleyFA, 1)
FRAC_MFT(vector, float, mworleyFB, 1)
FRAC_MFT(vector, float, cworleyFA, 1)
FRAC_MFT(vector, float, cworleyFB, 1)
FRAC_MFT(vector, float, alligator, 1)

// Definitions for: mfT, 2D
//------------------------------------------------------------------------------
FRAC_MFT( float, vector2, sine,      2)
FRAC_MFT( float, vector2, perlin,    2)
FRAC_MFT( float, vector2, pperlin,   2)
FRAC_MFT( float, vector2, simplex,   2)
FRAC_MFT( float, vector2, sparse,    2)
FRAC_MFT( float, vector2, flow,      2)
FRAC_MFT( float, vector2, pflow,     2)
FRAC_MFT( float, vector2, worleyFA,  2)
FRAC_MFT( float, vector2, worleyFB,  2)
FRAC_MFT( float, vector2, mworleyFA, 2)
FRAC_MFT( float, vector2, mworleyFB, 2)
FRAC_MFT( float, vector2, cworleyFA, 2)
FRAC_MFT( float, vector2, cworleyFB, 2)
FRAC_MFT( float, vector2, alligator, 2)
FRAC_MFT(vector, vector2, sine,      2)
FRAC_MFT(vector, vector2, perlin,    2)
FRAC_MFT(vector, vector2, pperlin,   2)
FRAC_MFT(vector, vector2, simplex,   2)
FRAC_MFT(vector, vector2, sparse,    2)
FRAC_MFT(vector, vector2, flow,      2)
FRAC_MFT(vector, vector2, pflow,     2)
FRAC_MFT(vector, vector2, worleyFA,  2)
FRAC_MFT(vector, vector2, worleyFB,  2)
FRAC_MFT(vector, vector2, mworleyFA, 2)
FRAC_MFT(vector, vector2, mworleyFB, 2)
FRAC_MFT(vector, vector2, cworleyFA, 2)
FRAC_MFT(vector, vector2, cworleyFB, 2)
FRAC_MFT(vector, vector2, alligator, 2)

// Definitions for: mfT, 3D
//------------------------------------------------------------------------------
FRAC_MFT( float, vector, sine,      3)
FRAC_MFT( float, vector, perlin,    3)
FRAC_MFT( float, vector, pperlin,   3)
FRAC_MFT( float, vector, simplex,   3)
FRAC_MFT( float, vector, sparse,    3)
FRAC_MFT( float, vector, flow,      3)
FRAC_MFT( float, vector, pflow,     3)
FRAC_MFT( float, vector, worleyFA,  3)
FRAC_MFT( float, vector, worleyFB,  3)
FRAC_MFT( float, vector, mworleyFA, 3)
FRAC_MFT( float, vector, mworleyFB, 3)
FRAC_MFT( float, vector, cworleyFA, 3)
FRAC_MFT( float, vector, cworleyFB, 3)
FRAC_MFT( float, vector, alligator, 3)
FRAC_MFT(vector, vector, sine,      3)
FRAC_MFT(vector, vector, perlin,    3)
FRAC_MFT(vector, vector, pperlin,   3)
FRAC_MFT(vector, vector, simplex,   3)
FRAC_MFT(vector, vector, sparse,    3)
FRAC_MFT(vector, vector, flow,      3)
FRAC_MFT(vector, vector, pflow,     3)
FRAC_MFT(vector, vector, worleyFA,  3)
FRAC_MFT(vector, vector, worleyFB,  3)
FRAC_MFT(vector, vector, mworleyFA, 3)
FRAC_MFT(vector, vector, mworleyFB, 3)
FRAC_MFT(vector, vector, cworleyFA, 3)
FRAC_MFT(vector, vector, cworleyFB, 3)
FRAC_MFT(vector, vector, alligator, 3)

// Definitions for: mfT, 4D
//------------------------------------------------------------------------------
FRAC_MFT( float, vector4, sine,      4)
FRAC_MFT( float, vector4, perlin,    4)
FRAC_MFT( float, vector4, pperlin,   4)
FRAC_MFT( float, vector4, simplex,   4)
FRAC_MFT( float, vector4, sparse,    4)
FRAC_MFT( float, vector4, flow,      4)
FRAC_MFT( float, vector4, pflow,     4)
FRAC_MFT( float, vector4, worleyFA,  4)
FRAC_MFT( float, vector4, worleyFB,  4)
FRAC_MFT( float, vector4, mworleyFA, 4)
FRAC_MFT( float, vector4, mworleyFB, 4)
FRAC_MFT( float, vector4, cworleyFA, 4)
FRAC_MFT( float, vector4, cworleyFB, 4)
FRAC_MFT( float, vector4, alligator, 4)
FRAC_MFT(vector, vector4, sine,      4)
FRAC_MFT(vector, vector4, perlin,    4)
FRAC_MFT(vector, vector4, pperlin,   4)
FRAC_MFT(vector, vector4, simplex,   4)
FRAC_MFT(vector, vector4, sparse,    4)
FRAC_MFT(vector, vector4, flow,      4)
FRAC_MFT(vector, vector4, pflow,     4)
FRAC_MFT(vector, vector4, worleyFA,  4)
FRAC_MFT(vector, vector4, worleyFB,  4)
FRAC_MFT(vector, vector4, mworleyFA, 4)
FRAC_MFT(vector, vector4, mworleyFB, 4)
FRAC_MFT(vector, vector4, cworleyFA, 4)
FRAC_MFT(vector, vector4, cworleyFB, 4)
FRAC_MFT(vector, vector4, alligator, 4)
#undef FRAC_MFT


// -----------------------------------------------
// Hybrid Multifractal Terrain (hmfT)
// Successive frequencies are scaled by the 
// function's local current and previous values
// -----------------------------------------------
#define FRAC_HMFT(RTYPE,PTYPE,NAME,DIM) \
RTYPE FNAME(hmfT,NAME) (PTYPE p; BPARMS(fsize,PTYPE); BXPARMS; \
                  FPARMS; FXPARMS ) \
{ \
   nsdata stats = SNAME##RTYPE(NAME,DIM); \
   nndata norm  = nn_##NAME; \
   string name  = stats.name; \
   float  oct   = x_oct = 0; \
   float  med   = x_avg = navg_correct(stats,BARGS(fsize)); \
   float  pvt   = navg_correct(npivot_##NAME(stats),stats,BARGS(fsize)); \
   RTYPE  base  = BNAME(NAME) (p,BARGS(fsize),BXARGS); \
   if(fsize>=1 || oct>maxoctaves || lacunarity*maxoctaves*gain==0) \
      return base; \
\
   PTYPE   pp    = p; \
   RTYPE   sbase = base - pvt; \
   float   lac   = abs(lacunarity); \
   float   g     = abs(gain) * min(1,lac); \
   float   fw    = fsize*(1+fold), fwold;\
   RTYPE   out   = sbase, chaos; \
   RTYPE   wsum  = 1, w = base*2; \
   RTYPE   wk    = lerp(10, fold*stats.symmetric ? 10 : norm.scale, \
                        smooth(1,norm.blend,lac)); \
\
   do { \
      w *= g; \
      oct += 1; \
      if(oct>maxoctaves) { \
         oct = maxoctaves; \
         w  *= bias(frac(maxoctaves),0.5-0.05); \
      } \
      fwold = fw; pp *= lac; fw *= lac; \
      if(fw>1) { \
         float k = (1-fwold) / (fw-fwold); \
         oct -= 1-k; \
         w   *= k; \
      } \
      chaos = BNAME(NAME) (pp,BARGS(fsize),BXARGS); \
      chaos -= pvt; \
      out  += w * chaos; \
      wsum += wk * w*w; \
      w *= 2*(chaos + pvt); \
\
   } while(oct<maxoctaves && fw<1.0); \
\
   x_oct = oct; \
   if(wsum!=0) out /= sqrt(wsum); \
   out = out*fit(g,0,0.1,1,natten(out+pvt,0.9,2)) + pvt; \
   out = clamp(out,(RTYPE)0,(RTYPE)1); /* safety - there should be no clipping */\
   return out; \
}


// Definitions for: hmfT, 1D
//------------------------------------------------------------------------------
FRAC_HMFT( float, float, sine,      1)
FRAC_HMFT( float, float, perlin,    1)
FRAC_HMFT( float, float, pperlin,   1)
FRAC_HMFT( float, float, simplex,   1)
FRAC_HMFT( float, float, sparse,    1)
FRAC_HMFT( float, float, flow,      1)
FRAC_HMFT( float, float, pflow,     1)
FRAC_HMFT( float, float, worleyFA,  1)
FRAC_HMFT( float, float, worleyFB,  1)
FRAC_HMFT( float, float, mworleyFA, 1)
FRAC_HMFT( float, float, mworleyFB, 1)
FRAC_HMFT( float, float, cworleyFA, 1)
FRAC_HMFT( float, float, cworleyFB, 1)
FRAC_HMFT( float, float, alligator, 1)
FRAC_HMFT(vector, float, sine,      1)
FRAC_HMFT(vector, float, perlin,    1)
FRAC_HMFT(vector, float, pperlin,   1)
FRAC_HMFT(vector, float, simplex,   1)
FRAC_HMFT(vector, float, sparse,    1)
FRAC_HMFT(vector, float, flow,      1)
FRAC_HMFT(vector, float, pflow,     1)
FRAC_HMFT(vector, float, worleyFA,  1)
FRAC_HMFT(vector, float, worleyFB,  1)
FRAC_HMFT(vector, float, mworleyFA, 1)
FRAC_HMFT(vector, float, mworleyFB, 1)
FRAC_HMFT(vector, float, cworleyFA, 1)
FRAC_HMFT(vector, float, cworleyFB, 1)
FRAC_HMFT(vector, float, alligator, 1)

// Definitions for: hmfT, 2D
//------------------------------------------------------------------------------
FRAC_HMFT( float, vector2, sine,      2)
FRAC_HMFT( float, vector2, perlin,    2)
FRAC_HMFT( float, vector2, pperlin,   2)
FRAC_HMFT( float, vector2, simplex,   2)
FRAC_HMFT( float, vector2, sparse,    2)
FRAC_HMFT( float, vector2, flow,      2)
FRAC_HMFT( float, vector2, pflow,     2)
FRAC_HMFT( float, vector2, worleyFA,  2)
FRAC_HMFT( float, vector2, worleyFB,  2)
FRAC_HMFT( float, vector2, mworleyFA, 2)
FRAC_HMFT( float, vector2, mworleyFB, 2)
FRAC_HMFT( float, vector2, cworleyFA, 2)
FRAC_HMFT( float, vector2, cworleyFB, 2)
FRAC_HMFT( float, vector2, alligator, 2)
FRAC_HMFT(vector, vector2, sine,      2)
FRAC_HMFT(vector, vector2, perlin,    2)
FRAC_HMFT(vector, vector2, pperlin,   2)
FRAC_HMFT(vector, vector2, simplex,   2)
FRAC_HMFT(vector, vector2, sparse,    2)
FRAC_HMFT(vector, vector2, flow,      2)
FRAC_HMFT(vector, vector2, pflow,     2)
FRAC_HMFT(vector, vector2, worleyFA,  2)
FRAC_HMFT(vector, vector2, worleyFB,  2)
FRAC_HMFT(vector, vector2, mworleyFA, 2)
FRAC_HMFT(vector, vector2, mworleyFB, 2)
FRAC_HMFT(vector, vector2, cworleyFA, 2)
FRAC_HMFT(vector, vector2, cworleyFB, 2)
FRAC_HMFT(vector, vector2, alligator, 2)

// Definitions for: hmfT, 3D
//------------------------------------------------------------------------------
FRAC_HMFT( float, vector, sine,      3)
FRAC_HMFT( float, vector, perlin,    3)
FRAC_HMFT( float, vector, pperlin,   3)
FRAC_HMFT( float, vector, simplex,   3)
FRAC_HMFT( float, vector, sparse,    3)
FRAC_HMFT( float, vector, flow,      3)
FRAC_HMFT( float, vector, pflow,     3)
FRAC_HMFT( float, vector, worleyFA,  3)
FRAC_HMFT( float, vector, worleyFB,  3)
FRAC_HMFT( float, vector, mworleyFA, 3)
FRAC_HMFT( float, vector, mworleyFB, 3)
FRAC_HMFT( float, vector, cworleyFA, 3)
FRAC_HMFT( float, vector, cworleyFB, 3)
FRAC_HMFT( float, vector, alligator, 3)
FRAC_HMFT(vector, vector, sine,      3)
FRAC_HMFT(vector, vector, perlin,    3)
FRAC_HMFT(vector, vector, pperlin,   3)
FRAC_HMFT(vector, vector, simplex,   3)
FRAC_HMFT(vector, vector, sparse,    3)
FRAC_HMFT(vector, vector, flow,      3)
FRAC_HMFT(vector, vector, pflow,     3)
FRAC_HMFT(vector, vector, worleyFA,  3)
FRAC_HMFT(vector, vector, worleyFB,  3)
FRAC_HMFT(vector, vector, mworleyFA, 3)
FRAC_HMFT(vector, vector, mworleyFB, 3)
FRAC_HMFT(vector, vector, cworleyFA, 3)
FRAC_HMFT(vector, vector, cworleyFB, 3)
FRAC_HMFT(vector, vector, alligator, 3)

// Definitions for: hmfT, 4D
//------------------------------------------------------------------------------
FRAC_HMFT( float, vector4, sine,      4)
FRAC_HMFT( float, vector4, perlin,    4)
FRAC_HMFT( float, vector4, pperlin,   4)
FRAC_HMFT( float, vector4, simplex,   4)
FRAC_HMFT( float, vector4, sparse,    4)
FRAC_HMFT( float, vector4, flow,      4)
FRAC_HMFT( float, vector4, pflow,     4)
FRAC_HMFT( float, vector4, worleyFA,  4)
FRAC_HMFT( float, vector4, worleyFB,  4)
FRAC_HMFT( float, vector4, mworleyFA, 4)
FRAC_HMFT( float, vector4, mworleyFB, 4)
FRAC_HMFT( float, vector4, cworleyFA, 4)
FRAC_HMFT( float, vector4, cworleyFB, 4)
FRAC_HMFT( float, vector4, alligator, 4)
FRAC_HMFT(vector, vector4, sine,      4)
FRAC_HMFT(vector, vector4, perlin,    4)
FRAC_HMFT(vector, vector4, pperlin,   4)
FRAC_HMFT(vector, vector4, simplex,   4)
FRAC_HMFT(vector, vector4, sparse,    4)
FRAC_HMFT(vector, vector4, flow,      4)
FRAC_HMFT(vector, vector4, pflow,     4)
FRAC_HMFT(vector, vector4, worleyFA,  4)
FRAC_HMFT(vector, vector4, worleyFB,  4)
FRAC_HMFT(vector, vector4, mworleyFA, 4)
FRAC_HMFT(vector, vector4, mworleyFB, 4)
FRAC_HMFT(vector, vector4, cworleyFA, 4)
FRAC_HMFT(vector, vector4, cworleyFB, 4)
FRAC_HMFT(vector, vector4, alligator, 4)
#undef FRAC_HMFT




//------------------------------------------------------------------------------
// Output Correction
//------------------------------------------------------------------------------
float noise_cc(
      float  chaos;      // noise value to correct
      int    inv;        // whether to take the complement
      int    dobias;     // whether to apply bias corr.
      float  ccbias;     // bias [0,1]
      int    dogain;     // whether to apply gain corr.
      float  ccgain;     // gain [0,1]     
      int    dorng;      // whether to apply range mapping
      float  rnglo;      // ouput range: min [-inf,rnghi]
      float  rnghi;      // ouput range: max [rnglo,+inf]
      float  amp;        // final amplitude
   )
{
   float out = chaos;
   if(dobias) out = bias(out,clamp(ccbias,1e-5,1-1e-5));
   if(dogain) out = gain(out,clamp(ccgain,1e-5,1-1e-5));
   if(inv) out = 1-out;
   if(dorng) out = fit(out,0,1,rnglo,rnghi);
   return out*amp;
}
vector noise_cc(
      vector chaos;      // noise value to correct
      int    inv;        // whether to take the complement
      int    dobias;     // whether to apply bias corr.
      float  ccbias;     // bias [0,1]
      int    dogain;     // whether to apply gain corr.
      float  ccgain;     // gain [0,1]     
      int    dorng;      // whether to apply range mapping
      float  rnglo;      // ouput range: min [-inf,rnghi]
      float  rnghi;      // ouput range: max [rnglo,+inf]
      float  amp;        // final amplitude
   )
{
   vector out = chaos;
   if(dobias) out = bias(out,clamp(ccbias,1e-5,1-1e-5));
   if(dogain) out = gain(out,clamp(ccgain,1e-5,1-1e-5));
   if(inv) out = VONE-out;
   if(dorng) out = fit(out,VZERO,VONE,(vector)rnglo,(vector)rnghi);
   return out*amp;
}

float noise_cc(
      float  chaos;      // noise value to correct
      int    inv;        // whether to take the complement
      int    dobias;     // whether to apply bias corr.
      vector ccbias;     // bias [0,1]
      int    dogain;     // whether to apply gain corr.
      vector ccgain;     // gain [0,1]     
      int    dorng;      // whether to apply range mapping
      vector rnglo;      // ouput range: min [-inf,rnghi]
      vector rnghi;      // ouput range: max [rnglo,+inf]
      vector amp;        // final amplitude
   )
{
   float out = chaos;
   if(dobias) out = bias(out,clamp(ccbias[0],1e-5,1-1e-5));
   if(dogain) out = gain(out,clamp(ccgain[0],1e-5,1-1e-5));
   if(inv) out = 1-out;
   if(dorng) out = fit(out,0,1,rnglo[0],rnghi[0]);
   return out*amp[0];
}
vector noise_cc(
      vector chaos;      // noise value to correct
      int    inv;        // whether to take the complement
      int    dobias;     // whether to apply bias corr.
      vector ccbias;     // bias [0,1]
      int    dogain;     // whether to apply gain corr.
      vector ccgain;     // gain [0,1]     
      int    dorng;      // whether to apply range mapping
      vector rnglo;      // ouput range: min [-inf,rnghi]
      vector rnghi;      // ouput range: max [rnglo,+inf]
      vector amp;        // final amplitude
   )
{
   vector out = chaos;
   if(dobias) out = bias(out,clamp(ccbias,(vector)1e-5,(vector)1-1e-5));
   if(dogain) out = gain(out,clamp(ccgain,(vector)1e-5,(vector)1-1e-5));
   if(inv) out = VONE-out;
   if(dorng) out = fit(out,VZERO,VONE,rnglo,rnghi);
   return out*amp;
}

// FOR BACKWARDS COMPATIBILITY WITH <= 14.0
float fnoise_cc(
      float  chaos;
      int    inv;
      int    dobias;
      float  ccbias;
      int    dogain;
      float  ccgain;
      int    dorng;
      float  rnglo;
      float  rnghi;
      float  amp;
   )
{
   return noise_cc(chaos, inv, dobias, ccbias, dogain, ccgain, dorng, rnglo, rnghi, amp);
}

vector vnoise_cc(
      float  chaos;
      int    inv;
      int    dobias;
      float  ccbias;
      int    dogain;
      float  ccgain;
      int    dorng;
      float  rnglo;
      float  rnghi;
      float  amp;
   )
{
   return noise_cc(chaos, inv, dobias, ccbias, dogain, ccgain, dorng, rnglo, rnghi, amp);
}

float fnoise_ccv(
      float  chaos;
      int    inv;
      int    dobias;
      vector  ccbias;
      int    dogain;
      vector  ccgain;
      int    dorng;
      vector  rnglo;
      vector  rnghi;
      vector  amp;
   )
{
   return noise_cc(chaos, inv, dobias, ccbias, dogain, ccgain, dorng, rnglo, rnghi, amp);
}

vector vnoise_ccv(
      float  chaos;
      int    inv;
      int    dobias;
      vector  ccbias;
      int    dogain;
      vector  ccgain;
      int    dorng;
      vector  rnglo;
      vector  rnghi;
      vector  amp;
   )
{
   return noise_cc(chaos, inv, dobias, ccbias, dogain, ccgain, dorng, rnglo, rnghi, amp);
}

#define UNIFIEDBASE(RTYPE, PTYPE) \
RTYPE unified_noise(string basis; PTYPE p; BPARMS(fw, PTYPE); BXPARMS) \
{ \
\
    if(basis == 'sine') \
        return noise_sine(p, BARGS(fw), BXARGS); \
    else if(basis == 'perlin') \
        return noise_perlin(p, BARGS(fw), BXARGS); \
    else if(basis == 'pperlin') \
        return noise_pperlin(p, BARGS(fw), BXARGS); \
    else if(basis == 'simplex') \
        return noise_simplex(p, BARGS(fw), BXARGS); \
    else if(basis == 'sparse') \
        return noise_sparse(p, BARGS(fw), BXARGS); \
    else if(basis == 'flow') \
        return noise_flow(p, BARGS(fw), BXARGS); \
    else if(basis == 'pflow') \
        return noise_pflow(p, BARGS(fw), BXARGS); \
    else if(basis == 'worleyFA') \
        return noise_worleyFA(p, BARGS(fw), BXARGS); \
    else if(basis == 'worleyFB') \
        return noise_worleyFB(p, BARGS(fw), BXARGS); \
    else if(basis == 'mworleyFA') \
        return noise_mworleyFA(p, BARGS(fw), BXARGS); \
    else if(basis == 'mworleyFB') \
        return noise_mworleyFB(p, BARGS(fw), BXARGS); \
    else if(basis == 'cworleyFA') \
        return noise_cworleyFA(p, BARGS(fw), BXARGS); \
    else if(basis == 'cworleyFB') \
        return noise_cworleyFB(p, BARGS(fw), BXARGS); \
    else if(basis == 'alligator') \
        return noise_alligator(p, BARGS(fw), BXARGS); \
\
    return 1.0; \
}

UNIFIEDBASE(float, float)
UNIFIEDBASE(float, vector2)
UNIFIEDBASE(float, vector)
UNIFIEDBASE(float, vector4)
UNIFIEDBASE(vector, float)
UNIFIEDBASE(vector, vector2)
UNIFIEDBASE(vector, vector)
UNIFIEDBASE(vector, vector4)

#define FRACBLOCK(BASIS) \
if(fractal == 'fBm') \
    return fBm_##BASIS(p, BARGS(fw), BXARGS, FARGS, FXARGS); \
else if(fractal == 'mfT') \
    return mfT_##BASIS(p, BARGS(fw), BXARGS, FARGS, FXARGS); \
else if(fractal == 'hmfT') \
    return hmfT_##BASIS(p, BARGS(fw), BXARGS, FARGS, FXARGS); \
else \
    return 1.0;

#define UNIFIEDFRACTAL(RTYPE, PTYPE) \
RTYPE unified_fractal_noise(string fractal, basis; PTYPE p; \
		    BPARMS(fw, PTYPE); BXPARMS; \
		    FPARMS; FXPARMS) \
{ \
\
    if(basis == 'sine') \
        FRACBLOCK(sine) \
    else if(basis == 'perlin') \
        FRACBLOCK(perlin) \
    else if(basis == 'pperlin') \
        FRACBLOCK(pperlin) \
    else if(basis == 'simplex') \
        FRACBLOCK(simplex) \
    else if(basis == 'sparse') \
        FRACBLOCK(sparse) \
    else if(basis == 'flow') \
        FRACBLOCK(flow) \
    else if(basis == 'pflow') \
        FRACBLOCK(pflow) \
    else if(basis == 'worleyFA') \
        FRACBLOCK(worleyFA) \
    else if(basis == 'worleyFB') \
        FRACBLOCK(worleyFB) \
    else if(basis == 'mworleyFA') \
        FRACBLOCK(mworleyFA) \
    else if(basis == 'mworleyFB') \
        FRACBLOCK(mworleyFB) \
    else if(basis == 'cworleyFA') \
        FRACBLOCK(cworleyFA) \
    else if(basis == 'cworleyFB') \
        FRACBLOCK(cworleyFB) \
    else if(basis == 'alligator') \
        FRACBLOCK(alligator) \
    return 1.0; \
}

UNIFIEDFRACTAL(float, float)
UNIFIEDFRACTAL(float, vector2)
UNIFIEDFRACTAL(float, vector)
UNIFIEDFRACTAL(float, vector4)
UNIFIEDFRACTAL(vector, float)
UNIFIEDFRACTAL(vector, vector2)
UNIFIEDFRACTAL(vector, vector)
UNIFIEDFRACTAL(vector, vector4)

// FOR BACKWARDS COMPATIBILITY WITH <= 14.0
// Keep functions with v/f prefix and dimension postfix. The new functions
// above don't need these because they use function overloading

#define BNAMEOLDfloat(name,dim)       CAT3(fnoise_,name,dim)
#define BNAMEOLDvector(name,dim)      CAT3(vnoise_,name,dim)

#define BARGS_OLD(FW)  inv,fold,accl,accg,FW,expon,_per, \
                   dolw,dogw,disp,dfreq,gflow,flow

#define DEF_BASE_OLD(RTYPE,PTYPE,NAME,DIM) \
RTYPE BNAMEOLD##RTYPE(NAME,DIM) (vector4 p; BPARMS(fsize, vector4); BXPARMS) \
{ \
    PTYPE _p = (PTYPE) p; \
    PTYPE _per = (PTYPE) per; \
    RTYPE res = BNAME(NAME) (_p, BARGS_OLD(fsize), BXARGS); \
    return (RTYPE) res; \
}

// Old names & signatures for backwards compatibility with <= 14.0

DEF_BASE_OLD(float ,float,  sine,1)
DEF_BASE_OLD(float ,vector2,sine,2)
DEF_BASE_OLD(float ,vector, sine,3)
DEF_BASE_OLD(float ,vector4,sine,4)
DEF_BASE_OLD(vector,float,  sine,1)
DEF_BASE_OLD(vector,vector2,sine,2)
DEF_BASE_OLD(vector,vector, sine,3)
DEF_BASE_OLD(vector,vector4,sine,4)

DEF_BASE_OLD(float ,float,  perlin,1)
DEF_BASE_OLD(float ,vector2,perlin,2)
DEF_BASE_OLD(float ,vector, perlin,3)
DEF_BASE_OLD(float ,vector4,perlin,4)
DEF_BASE_OLD(vector,float,  perlin,1)
DEF_BASE_OLD(vector,vector2,perlin,2)
DEF_BASE_OLD(vector,vector, perlin,3)
DEF_BASE_OLD(vector,vector4,perlin,4)

DEF_BASE_OLD(float ,float,  pperlin,1)
DEF_BASE_OLD(float ,vector2,pperlin,2)
DEF_BASE_OLD(float ,vector, pperlin,3)
DEF_BASE_OLD(float ,vector4,pperlin,4)
DEF_BASE_OLD(vector,float,  pperlin,1)
DEF_BASE_OLD(vector,vector2,pperlin,2)
DEF_BASE_OLD(vector,vector, pperlin,3)
DEF_BASE_OLD(vector,vector4,pperlin,4)

DEF_BASE_OLD(float ,float,  simplex,1)
DEF_BASE_OLD(float ,vector2,simplex,2)
DEF_BASE_OLD(float ,vector, simplex,3)
DEF_BASE_OLD(float ,vector4,simplex,4)
DEF_BASE_OLD(vector,float,  simplex,1)
DEF_BASE_OLD(vector,vector2,simplex,2)
DEF_BASE_OLD(vector,vector, simplex,3)
DEF_BASE_OLD(vector,vector4,simplex,4)

DEF_BASE_OLD(float ,float,  sparse,1)
DEF_BASE_OLD(float ,vector2,sparse,2)
DEF_BASE_OLD(float ,vector, sparse,3)
DEF_BASE_OLD(float ,vector4,sparse,4)
DEF_BASE_OLD(vector,float,  sparse,1)
DEF_BASE_OLD(vector,vector2,sparse,2)
DEF_BASE_OLD(vector,vector, sparse,3)
DEF_BASE_OLD(vector,vector4,sparse,4)

DEF_BASE_OLD(float ,float,  flow,1)
DEF_BASE_OLD(float ,vector2,flow,2)
DEF_BASE_OLD(float ,vector, flow,3)
DEF_BASE_OLD(float ,vector4,flow,4)
DEF_BASE_OLD(vector,float,  flow,1)
DEF_BASE_OLD(vector,vector2,flow,2)
DEF_BASE_OLD(vector,vector, flow,3)
DEF_BASE_OLD(vector,vector4,flow,4)

DEF_BASE_OLD(float ,float,  pflow,1)
DEF_BASE_OLD(float ,vector2,pflow,2)
DEF_BASE_OLD(float ,vector, pflow,3)
DEF_BASE_OLD(float ,vector4,pflow,4)
DEF_BASE_OLD(vector,float,  pflow,1)
DEF_BASE_OLD(vector,vector2,pflow,2)
DEF_BASE_OLD(vector,vector, pflow,3)
DEF_BASE_OLD(vector,vector4,pflow,4)

DEF_BASE_OLD(float ,float,  worleyFA,1)
DEF_BASE_OLD(float ,vector2,worleyFA,2)
DEF_BASE_OLD(float ,vector, worleyFA,3)
DEF_BASE_OLD(float ,vector4,worleyFA,4)
DEF_BASE_OLD(vector,float,  worleyFA,1)
DEF_BASE_OLD(vector,vector2,worleyFA,2)
DEF_BASE_OLD(vector,vector, worleyFA,3)
DEF_BASE_OLD(vector,vector4,worleyFA,4)

DEF_BASE_OLD(float ,float,  worleyFB,1)
DEF_BASE_OLD(float ,vector2,worleyFB,2)
DEF_BASE_OLD(float ,vector, worleyFB,3)
DEF_BASE_OLD(float ,vector4,worleyFB,4)
DEF_BASE_OLD(vector,float,  worleyFB,1)
DEF_BASE_OLD(vector,vector2,worleyFB,2)
DEF_BASE_OLD(vector,vector, worleyFB,3)
DEF_BASE_OLD(vector,vector4,worleyFB,4)

DEF_BASE_OLD(float ,float,  mworleyFA,1)
DEF_BASE_OLD(float ,vector2,mworleyFA,2)
DEF_BASE_OLD(float ,vector, mworleyFA,3)
DEF_BASE_OLD(float ,vector4,mworleyFA,4)
DEF_BASE_OLD(vector,float,  mworleyFA,1)
DEF_BASE_OLD(vector,vector2,mworleyFA,2)
DEF_BASE_OLD(vector,vector, mworleyFA,3)
DEF_BASE_OLD(vector,vector4,mworleyFA,4)

DEF_BASE_OLD(float ,float,  mworleyFB,1)
DEF_BASE_OLD(float ,vector2,mworleyFB,2)
DEF_BASE_OLD(float ,vector, mworleyFB,3)
DEF_BASE_OLD(float ,vector4,mworleyFB,4)
DEF_BASE_OLD(vector,float,  mworleyFB,1)
DEF_BASE_OLD(vector,vector2,mworleyFB,2)
DEF_BASE_OLD(vector,vector, mworleyFB,3)
DEF_BASE_OLD(vector,vector4,mworleyFB,4)

DEF_BASE_OLD(float ,float,  cworleyFA,1)
DEF_BASE_OLD(float ,vector2,cworleyFA,2)
DEF_BASE_OLD(float ,vector, cworleyFA,3)
DEF_BASE_OLD(float ,vector4,cworleyFA,4)
DEF_BASE_OLD(vector,float,  cworleyFA,1)
DEF_BASE_OLD(vector,vector2,cworleyFA,2)
DEF_BASE_OLD(vector,vector, cworleyFA,3)
DEF_BASE_OLD(vector,vector4,cworleyFA,4)

DEF_BASE_OLD(float ,float,  cworleyFB,1)
DEF_BASE_OLD(float ,vector2,cworleyFB,2)
DEF_BASE_OLD(float ,vector, cworleyFB,3)
DEF_BASE_OLD(float ,vector4,cworleyFB,4)
DEF_BASE_OLD(vector,float,  cworleyFB,1)
DEF_BASE_OLD(vector,vector2,cworleyFB,2)
DEF_BASE_OLD(vector,vector, cworleyFB,3)
DEF_BASE_OLD(vector,vector4,cworleyFB,4)

DEF_BASE_OLD(float ,float,  alligator,1)
DEF_BASE_OLD(float ,vector2,alligator,2)
DEF_BASE_OLD(float ,vector, alligator,3)
DEF_BASE_OLD(float ,vector4,alligator,4)
DEF_BASE_OLD(vector,float,  alligator,1)
DEF_BASE_OLD(vector,vector2,alligator,2)
DEF_BASE_OLD(vector,vector, alligator,3)
DEF_BASE_OLD(vector,vector4,alligator,4)
#undef DEF_BASE_OLD

// -------------------------------------------------------------------------------
// Fractal Constructions
// Naming: <frac_type>_<v|f>_<base_name>_<dim>
//     eg: fBm_fperlin1() or fBm_vperlin1()
// -------------------------------------------------------------------------------
#define FNAMEOLDfloat(pfx,name,dim)  CAT4(pfx##_,f,name,dim)
#define FNAMEOLDvector(pfx,name,dim) CAT4(pfx##_,v,name,dim)


#define BARGS_OLD(FW)  inv,fold,accl,accg,FW,expon,_per, \
                   dolw,dogw,disp,dfreq,gflow,flow


#define FRAC_OLD(RTYPE,FRACTYPE,PTYPE,NAME,DIM) \
RTYPE FNAMEOLD##RTYPE(FRACTYPE,NAME,DIM) (vector4 p; BPARMS(fsize, vector4); BXPARMS; \
                  FPARMS; FXPARMS ) \
{ \
    PTYPE _p = (PTYPE) p; \
    PTYPE _per = (PTYPE) per; \
    RTYPE res = FNAME(FRACTYPE, NAME) (_p, BARGS_OLD(fsize), BXARGS, FARGS, FXARGS); \
    return (RTYPE) res; \
}

// Old names & signatures for backwards compatibility with <= 14.0

FRAC_OLD( float, fBm, float, sine,      1)
FRAC_OLD( float, fBm, float, perlin,    1)
FRAC_OLD( float, fBm, float, pperlin,   1)
FRAC_OLD( float, fBm, float, simplex,   1)
FRAC_OLD( float, fBm, float, sparse,    1)
FRAC_OLD( float, fBm, float, flow,      1)
FRAC_OLD( float, fBm, float, pflow,     1)
FRAC_OLD( float, fBm, float, worleyFA,  1)
FRAC_OLD( float, fBm, float, worleyFB,  1)
FRAC_OLD( float, fBm, float, mworleyFA, 1)
FRAC_OLD( float, fBm, float, mworleyFB, 1)
FRAC_OLD( float, fBm, float, cworleyFA, 1)
FRAC_OLD( float, fBm, float, cworleyFB, 1)
FRAC_OLD( float, fBm, float, alligator, 1)
FRAC_OLD(vector, fBm, float, sine,      1)
FRAC_OLD(vector, fBm, float, perlin,    1)
FRAC_OLD(vector, fBm, float, pperlin,   1)
FRAC_OLD(vector, fBm, float, simplex,   1)
FRAC_OLD(vector, fBm, float, sparse,    1)
FRAC_OLD(vector, fBm, float, flow,      1)
FRAC_OLD(vector, fBm, float, pflow,     1)
FRAC_OLD(vector, fBm, float, worleyFA,  1)
FRAC_OLD(vector, fBm, float, worleyFB,  1)
FRAC_OLD(vector, fBm, float, mworleyFA, 1)
FRAC_OLD(vector, fBm, float, mworleyFB, 1)
FRAC_OLD(vector, fBm, float, cworleyFA, 1)
FRAC_OLD(vector, fBm, float, cworleyFB, 1)
FRAC_OLD(vector, fBm, float, alligator, 1)

FRAC_OLD( float, fBm, vector2, perlin,    2)
FRAC_OLD( float, fBm, vector2, pperlin,   2)
FRAC_OLD( float, fBm, vector2, simplex,   2)
FRAC_OLD( float, fBm, vector2, sparse,    2)
FRAC_OLD( float, fBm, vector2, flow,      2)
FRAC_OLD( float, fBm, vector2, pflow,     2)
FRAC_OLD( float, fBm, vector2, worleyFA,  2)
FRAC_OLD( float, fBm, vector2, worleyFB,  2)
FRAC_OLD( float, fBm, vector2, mworleyFA, 2)
FRAC_OLD( float, fBm, vector2, mworleyFB, 2)
FRAC_OLD( float, fBm, vector2, cworleyFA, 2)
FRAC_OLD( float, fBm, vector2, cworleyFB, 2)
FRAC_OLD( float, fBm, vector2, alligator, 2)
FRAC_OLD(vector, fBm, vector2, sine,      2)
FRAC_OLD(vector, fBm, vector2, perlin,    2)
FRAC_OLD(vector, fBm, vector2, pperlin,   2)
FRAC_OLD(vector, fBm, vector2, simplex,   2)
FRAC_OLD(vector, fBm, vector2, sparse,    2)
FRAC_OLD(vector, fBm, vector2, flow,      2)
FRAC_OLD(vector, fBm, vector2, pflow,     2)
FRAC_OLD(vector, fBm, vector2, worleyFA,  2)
FRAC_OLD(vector, fBm, vector2, worleyFB,  2)
FRAC_OLD(vector, fBm, vector2, mworleyFA, 2)
FRAC_OLD(vector, fBm, vector2, mworleyFB, 2)
FRAC_OLD(vector, fBm, vector2, cworleyFA, 2)
FRAC_OLD(vector, fBm, vector2, cworleyFB, 2)
FRAC_OLD(vector, fBm, vector2, alligator, 2)

FRAC_OLD( float, fBm, vector, sine,      3)
FRAC_OLD( float, fBm, vector, perlin,    3)
FRAC_OLD( float, fBm, vector, pperlin,   3)
FRAC_OLD( float, fBm, vector, simplex,   3)
FRAC_OLD( float, fBm, vector, sparse,    3)
FRAC_OLD( float, fBm, vector, flow,      3)
FRAC_OLD( float, fBm, vector, pflow,     3)
FRAC_OLD( float, fBm, vector, worleyFA,  3)
FRAC_OLD( float, fBm, vector, worleyFB,  3)
FRAC_OLD( float, fBm, vector, mworleyFA, 3)
FRAC_OLD( float, fBm, vector, mworleyFB, 3)
FRAC_OLD( float, fBm, vector, cworleyFA, 3)
FRAC_OLD( float, fBm, vector, cworleyFB, 3)
FRAC_OLD( float, fBm, vector, alligator, 3)
FRAC_OLD(vector, fBm, vector, sine,      3)
FRAC_OLD(vector, fBm, vector, perlin,    3)
FRAC_OLD(vector, fBm, vector, pperlin,   3)
FRAC_OLD(vector, fBm, vector, simplex,   3)
FRAC_OLD(vector, fBm, vector, sparse,    3)
FRAC_OLD(vector, fBm, vector, flow,      3)
FRAC_OLD(vector, fBm, vector, pflow,     3)
FRAC_OLD(vector, fBm, vector, worleyFA,  3)
FRAC_OLD(vector, fBm, vector, worleyFB,  3)
FRAC_OLD(vector, fBm, vector, mworleyFA, 3)
FRAC_OLD(vector, fBm, vector, mworleyFB, 3)
FRAC_OLD(vector, fBm, vector, cworleyFA, 3)
FRAC_OLD(vector, fBm, vector, cworleyFB, 3)
FRAC_OLD(vector, fBm, vector, alligator, 3)

FRAC_OLD( float, fBm, vector4, sine,      4)
FRAC_OLD( float, fBm, vector4, perlin,    4)
FRAC_OLD( float, fBm, vector4, pperlin,   4)
FRAC_OLD( float, fBm, vector4, simplex,   4)
FRAC_OLD( float, fBm, vector4, sparse,    4)
FRAC_OLD( float, fBm, vector4, flow,      4)
FRAC_OLD( float, fBm, vector4, pflow,     4)
FRAC_OLD( float, fBm, vector4, worleyFA,  4)
FRAC_OLD( float, fBm, vector4, worleyFB,  4)
FRAC_OLD( float, fBm, vector4, mworleyFA, 4)
FRAC_OLD( float, fBm, vector4, mworleyFB, 4)
FRAC_OLD( float, fBm, vector4, cworleyFA, 4)
FRAC_OLD( float, fBm, vector4, cworleyFB, 4)
FRAC_OLD( float, fBm, vector4, alligator, 4)
FRAC_OLD(vector, fBm, vector4, sine,      4)
FRAC_OLD(vector, fBm, vector4, perlin,    4)
FRAC_OLD(vector, fBm, vector4, pperlin,   4)
FRAC_OLD(vector, fBm, vector4, simplex,   4)
FRAC_OLD(vector, fBm, vector4, sparse,    4)
FRAC_OLD(vector, fBm, vector4, flow,      4)
FRAC_OLD(vector, fBm, vector4, pflow,     4)
FRAC_OLD(vector, fBm, vector4, worleyFA,  4)
FRAC_OLD(vector, fBm, vector4, worleyFB,  4)
FRAC_OLD(vector, fBm, vector4, mworleyFA, 4)
FRAC_OLD(vector, fBm, vector4, mworleyFB, 4)
FRAC_OLD(vector, fBm, vector4, cworleyFA, 4)
FRAC_OLD(vector, fBm, vector4, cworleyFB, 4)
FRAC_OLD(vector, fBm, vector4, alligator, 4)

FRAC_OLD( float, mfT, float, sine,      1)
FRAC_OLD( float, mfT, float, perlin,    1)
FRAC_OLD( float, mfT, float, pperlin,   1)
FRAC_OLD( float, mfT, float, simplex,   1)
FRAC_OLD( float, mfT, float, sparse,    1)
FRAC_OLD( float, mfT, float, flow,      1)
FRAC_OLD( float, mfT, float, pflow,     1)
FRAC_OLD( float, mfT, float, worleyFA,  1)
FRAC_OLD( float, mfT, float, worleyFB,  1)
FRAC_OLD( float, mfT, float, mworleyFA, 1)
FRAC_OLD( float, mfT, float, mworleyFB, 1)
FRAC_OLD( float, mfT, float, cworleyFA, 1)
FRAC_OLD( float, mfT, float, cworleyFB, 1)
FRAC_OLD( float, mfT, float, alligator, 1)
FRAC_OLD(vector, mfT, float, sine,      1)
FRAC_OLD(vector, mfT, float, perlin,    1)
FRAC_OLD(vector, mfT, float, pperlin,   1)
FRAC_OLD(vector, mfT, float, simplex,   1)
FRAC_OLD(vector, mfT, float, sparse,    1)
FRAC_OLD(vector, mfT, float, flow,      1)
FRAC_OLD(vector, mfT, float, pflow,     1)
FRAC_OLD(vector, mfT, float, worleyFA,  1)
FRAC_OLD(vector, mfT, float, worleyFB,  1)
FRAC_OLD(vector, mfT, float, mworleyFA, 1)
FRAC_OLD(vector, mfT, float, mworleyFB, 1)
FRAC_OLD(vector, mfT, float, cworleyFA, 1)
FRAC_OLD(vector, mfT, float, cworleyFB, 1)
FRAC_OLD(vector, mfT, float, alligator, 1)

FRAC_OLD( float, mfT, vector2, perlin,    2)
FRAC_OLD( float, mfT, vector2, pperlin,   2)
FRAC_OLD( float, mfT, vector2, simplex,   2)
FRAC_OLD( float, mfT, vector2, sparse,    2)
FRAC_OLD( float, mfT, vector2, flow,      2)
FRAC_OLD( float, mfT, vector2, pflow,     2)
FRAC_OLD( float, mfT, vector2, worleyFA,  2)
FRAC_OLD( float, mfT, vector2, worleyFB,  2)
FRAC_OLD( float, mfT, vector2, mworleyFA, 2)
FRAC_OLD( float, mfT, vector2, mworleyFB, 2)
FRAC_OLD( float, mfT, vector2, cworleyFA, 2)
FRAC_OLD( float, mfT, vector2, cworleyFB, 2)
FRAC_OLD( float, mfT, vector2, alligator, 2)
FRAC_OLD(vector, mfT, vector2, sine,      2)
FRAC_OLD(vector, mfT, vector2, perlin,    2)
FRAC_OLD(vector, mfT, vector2, pperlin,   2)
FRAC_OLD(vector, mfT, vector2, simplex,   2)
FRAC_OLD(vector, mfT, vector2, sparse,    2)
FRAC_OLD(vector, mfT, vector2, flow,      2)
FRAC_OLD(vector, mfT, vector2, pflow,     2)
FRAC_OLD(vector, mfT, vector2, worleyFA,  2)
FRAC_OLD(vector, mfT, vector2, worleyFB,  2)
FRAC_OLD(vector, mfT, vector2, mworleyFA, 2)
FRAC_OLD(vector, mfT, vector2, mworleyFB, 2)
FRAC_OLD(vector, mfT, vector2, cworleyFA, 2)
FRAC_OLD(vector, mfT, vector2, cworleyFB, 2)
FRAC_OLD(vector, mfT, vector2, alligator, 2)

FRAC_OLD( float, mfT, vector, sine,      3)
FRAC_OLD( float, mfT, vector, perlin,    3)
FRAC_OLD( float, mfT, vector, pperlin,   3)
FRAC_OLD( float, mfT, vector, simplex,   3)
FRAC_OLD( float, mfT, vector, sparse,    3)
FRAC_OLD( float, mfT, vector, flow,      3)
FRAC_OLD( float, mfT, vector, pflow,     3)
FRAC_OLD( float, mfT, vector, worleyFA,  3)
FRAC_OLD( float, mfT, vector, worleyFB,  3)
FRAC_OLD( float, mfT, vector, mworleyFA, 3)
FRAC_OLD( float, mfT, vector, mworleyFB, 3)
FRAC_OLD( float, mfT, vector, cworleyFA, 3)
FRAC_OLD( float, mfT, vector, cworleyFB, 3)
FRAC_OLD( float, mfT, vector, alligator, 3)
FRAC_OLD(vector, mfT, vector, sine,      3)
FRAC_OLD(vector, mfT, vector, perlin,    3)
FRAC_OLD(vector, mfT, vector, pperlin,   3)
FRAC_OLD(vector, mfT, vector, simplex,   3)
FRAC_OLD(vector, mfT, vector, sparse,    3)
FRAC_OLD(vector, mfT, vector, flow,      3)
FRAC_OLD(vector, mfT, vector, pflow,     3)
FRAC_OLD(vector, mfT, vector, worleyFA,  3)
FRAC_OLD(vector, mfT, vector, worleyFB,  3)
FRAC_OLD(vector, mfT, vector, mworleyFA, 3)
FRAC_OLD(vector, mfT, vector, mworleyFB, 3)
FRAC_OLD(vector, mfT, vector, cworleyFA, 3)
FRAC_OLD(vector, mfT, vector, cworleyFB, 3)
FRAC_OLD(vector, mfT, vector, alligator, 3)

FRAC_OLD( float, mfT, vector4, sine,      4)
FRAC_OLD( float, mfT, vector4, perlin,    4)
FRAC_OLD( float, mfT, vector4, pperlin,   4)
FRAC_OLD( float, mfT, vector4, simplex,   4)
FRAC_OLD( float, mfT, vector4, sparse,    4)
FRAC_OLD( float, mfT, vector4, flow,      4)
FRAC_OLD( float, mfT, vector4, pflow,     4)
FRAC_OLD( float, mfT, vector4, worleyFA,  4)
FRAC_OLD( float, mfT, vector4, worleyFB,  4)
FRAC_OLD( float, mfT, vector4, mworleyFA, 4)
FRAC_OLD( float, mfT, vector4, mworleyFB, 4)
FRAC_OLD( float, mfT, vector4, cworleyFA, 4)
FRAC_OLD( float, mfT, vector4, cworleyFB, 4)
FRAC_OLD( float, mfT, vector4, alligator, 4)
FRAC_OLD(vector, mfT, vector4, sine,      4)
FRAC_OLD(vector, mfT, vector4, perlin,    4)
FRAC_OLD(vector, mfT, vector4, pperlin,   4)
FRAC_OLD(vector, mfT, vector4, simplex,   4)
FRAC_OLD(vector, mfT, vector4, sparse,    4)
FRAC_OLD(vector, mfT, vector4, flow,      4)
FRAC_OLD(vector, mfT, vector4, pflow,     4)
FRAC_OLD(vector, mfT, vector4, worleyFA,  4)
FRAC_OLD(vector, mfT, vector4, worleyFB,  4)
FRAC_OLD(vector, mfT, vector4, mworleyFA, 4)
FRAC_OLD(vector, mfT, vector4, mworleyFB, 4)
FRAC_OLD(vector, mfT, vector4, cworleyFA, 4)
FRAC_OLD(vector, mfT, vector4, cworleyFB, 4)
FRAC_OLD(vector, mfT, vector4, alligator, 4)

FRAC_OLD( float, hmfT, float, sine,      1)
FRAC_OLD( float, hmfT, float, perlin,    1)
FRAC_OLD( float, hmfT, float, pperlin,   1)
FRAC_OLD( float, hmfT, float, simplex,   1)
FRAC_OLD( float, hmfT, float, sparse,    1)
FRAC_OLD( float, hmfT, float, flow,      1)
FRAC_OLD( float, hmfT, float, pflow,     1)
FRAC_OLD( float, hmfT, float, worleyFA,  1)
FRAC_OLD( float, hmfT, float, worleyFB,  1)
FRAC_OLD( float, hmfT, float, mworleyFA, 1)
FRAC_OLD( float, hmfT, float, mworleyFB, 1)
FRAC_OLD( float, hmfT, float, cworleyFA, 1)
FRAC_OLD( float, hmfT, float, cworleyFB, 1)
FRAC_OLD( float, hmfT, float, alligator, 1)
FRAC_OLD(vector, hmfT, float, sine,      1)
FRAC_OLD(vector, hmfT, float, perlin,    1)
FRAC_OLD(vector, hmfT, float, pperlin,   1)
FRAC_OLD(vector, hmfT, float, simplex,   1)
FRAC_OLD(vector, hmfT, float, sparse,    1)
FRAC_OLD(vector, hmfT, float, flow,      1)
FRAC_OLD(vector, hmfT, float, pflow,     1)
FRAC_OLD(vector, hmfT, float, worleyFA,  1)
FRAC_OLD(vector, hmfT, float, worleyFB,  1)
FRAC_OLD(vector, hmfT, float, mworleyFA, 1)
FRAC_OLD(vector, hmfT, float, mworleyFB, 1)
FRAC_OLD(vector, hmfT, float, cworleyFA, 1)
FRAC_OLD(vector, hmfT, float, cworleyFB, 1)
FRAC_OLD(vector, hmfT, float, alligator, 1)

FRAC_OLD( float, hmfT, vector2, perlin,    2)
FRAC_OLD( float, hmfT, vector2, pperlin,   2)
FRAC_OLD( float, hmfT, vector2, simplex,   2)
FRAC_OLD( float, hmfT, vector2, sparse,    2)
FRAC_OLD( float, hmfT, vector2, flow,      2)
FRAC_OLD( float, hmfT, vector2, pflow,     2)
FRAC_OLD( float, hmfT, vector2, worleyFA,  2)
FRAC_OLD( float, hmfT, vector2, worleyFB,  2)
FRAC_OLD( float, hmfT, vector2, mworleyFA, 2)
FRAC_OLD( float, hmfT, vector2, mworleyFB, 2)
FRAC_OLD( float, hmfT, vector2, cworleyFA, 2)
FRAC_OLD( float, hmfT, vector2, cworleyFB, 2)
FRAC_OLD( float, hmfT, vector2, alligator, 2)
FRAC_OLD(vector, hmfT, vector2, sine,      2)
FRAC_OLD(vector, hmfT, vector2, perlin,    2)
FRAC_OLD(vector, hmfT, vector2, pperlin,   2)
FRAC_OLD(vector, hmfT, vector2, simplex,   2)
FRAC_OLD(vector, hmfT, vector2, sparse,    2)
FRAC_OLD(vector, hmfT, vector2, flow,      2)
FRAC_OLD(vector, hmfT, vector2, pflow,     2)
FRAC_OLD(vector, hmfT, vector2, worleyFA,  2)
FRAC_OLD(vector, hmfT, vector2, worleyFB,  2)
FRAC_OLD(vector, hmfT, vector2, mworleyFA, 2)
FRAC_OLD(vector, hmfT, vector2, mworleyFB, 2)
FRAC_OLD(vector, hmfT, vector2, cworleyFA, 2)
FRAC_OLD(vector, hmfT, vector2, cworleyFB, 2)
FRAC_OLD(vector, hmfT, vector2, alligator, 2)

FRAC_OLD( float, hmfT, vector, sine,      3)
FRAC_OLD( float, hmfT, vector, perlin,    3)
FRAC_OLD( float, hmfT, vector, pperlin,   3)
FRAC_OLD( float, hmfT, vector, simplex,   3)
FRAC_OLD( float, hmfT, vector, sparse,    3)
FRAC_OLD( float, hmfT, vector, flow,      3)
FRAC_OLD( float, hmfT, vector, pflow,     3)
FRAC_OLD( float, hmfT, vector, worleyFA,  3)
FRAC_OLD( float, hmfT, vector, worleyFB,  3)
FRAC_OLD( float, hmfT, vector, mworleyFA, 3)
FRAC_OLD( float, hmfT, vector, mworleyFB, 3)
FRAC_OLD( float, hmfT, vector, cworleyFA, 3)
FRAC_OLD( float, hmfT, vector, cworleyFB, 3)
FRAC_OLD( float, hmfT, vector, alligator, 3)
FRAC_OLD(vector, hmfT, vector, sine,      3)
FRAC_OLD(vector, hmfT, vector, perlin,    3)
FRAC_OLD(vector, hmfT, vector, pperlin,   3)
FRAC_OLD(vector, hmfT, vector, simplex,   3)
FRAC_OLD(vector, hmfT, vector, sparse,    3)
FRAC_OLD(vector, hmfT, vector, flow,      3)
FRAC_OLD(vector, hmfT, vector, pflow,     3)
FRAC_OLD(vector, hmfT, vector, worleyFA,  3)
FRAC_OLD(vector, hmfT, vector, worleyFB,  3)
FRAC_OLD(vector, hmfT, vector, mworleyFA, 3)
FRAC_OLD(vector, hmfT, vector, mworleyFB, 3)
FRAC_OLD(vector, hmfT, vector, cworleyFA, 3)
FRAC_OLD(vector, hmfT, vector, cworleyFB, 3)
FRAC_OLD(vector, hmfT, vector, alligator, 3)

FRAC_OLD( float, hmfT, vector4, sine,      4)
FRAC_OLD( float, hmfT, vector4, perlin,    4)
FRAC_OLD( float, hmfT, vector4, pperlin,   4)
FRAC_OLD( float, hmfT, vector4, simplex,   4)
FRAC_OLD( float, hmfT, vector4, sparse,    4)
FRAC_OLD( float, hmfT, vector4, flow,      4)
FRAC_OLD( float, hmfT, vector4, pflow,     4)
FRAC_OLD( float, hmfT, vector4, worleyFA,  4)
FRAC_OLD( float, hmfT, vector4, worleyFB,  4)
FRAC_OLD( float, hmfT, vector4, mworleyFA, 4)
FRAC_OLD( float, hmfT, vector4, mworleyFB, 4)
FRAC_OLD( float, hmfT, vector4, cworleyFA, 4)
FRAC_OLD( float, hmfT, vector4, cworleyFB, 4)
FRAC_OLD( float, hmfT, vector4, alligator, 4)
FRAC_OLD(vector, hmfT, vector4, sine,      4)
FRAC_OLD(vector, hmfT, vector4, perlin,    4)
FRAC_OLD(vector, hmfT, vector4, pperlin,   4)
FRAC_OLD(vector, hmfT, vector4, simplex,   4)
FRAC_OLD(vector, hmfT, vector4, sparse,    4)
FRAC_OLD(vector, hmfT, vector4, flow,      4)
FRAC_OLD(vector, hmfT, vector4, pflow,     4)
FRAC_OLD(vector, hmfT, vector4, worleyFA,  4)
FRAC_OLD(vector, hmfT, vector4, worleyFB,  4)
FRAC_OLD(vector, hmfT, vector4, mworleyFA, 4)
FRAC_OLD(vector, hmfT, vector4, mworleyFB, 4)
FRAC_OLD(vector, hmfT, vector4, cworleyFA, 4)
FRAC_OLD(vector, hmfT, vector4, cworleyFB, 4)
FRAC_OLD(vector, hmfT, vector4, alligator, 4)

#undef FRAC_OLD


#endif // End pyro_noise2_h
