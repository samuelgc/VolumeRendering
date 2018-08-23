#include "Blackbody.h"
#include "Math.h"


double planckpower(double T) {
   return PLANCK_S * T*T*T*T;
}

void xyztochroma(double xyz[]) {
   double d = xyz[0] + xyz[1] + xyz[2];
   if(d == 0) return;
   divide(xyz, d);
}

double fit(double val, double min, double max) {
    val -= min;
    return val / (max - min);
}

void bbspline(double T, double xyz[]) {
   double  Tfirst = 1667;
   double  Tlast  = 25000;

   if(T<Tfirst) {
       double fit_val = fit(T, 0, Tfirst);
       xyz[0] = (0.564639 - 0.73469) * fit_val + 0.73469;
       xyz[1] = (0.402887 - 0.26531) * fit_val + 0.26531;
       xyz[2] = (0.0324743 - 3.92832e-21) * fit_val + 3.92832e-21;
   } else if(T>Tlast) {
       double fit_val = fit(T, Tlast, 1e6);
       xyz[0] = (0.238815 - 0.252473) * fit_val + 0.252473;
       xyz[1] = (0.239092 - 0.252255) * fit_val + 0.252255;
       xyz[2] = (0.522093 - 0.495272) * fit_val + 0.495272;
   } else {
      double x, y;
      double c1 = 1e9 / (T * T * T);
      double c2 = 1e6 / (T * T);
      double c3 = 1e3 / T; 
      if(T >= 1667 && T <= 4000) { 
         x = -0.2661239 * c1 - 0.2343580 * c2 + 0.8776956 * c3 + 0.179910; 
      } else if(T > 4000 && T <= 25000) { 
         x = -3.0258469 * c1 + 2.1070379 * c2 + 0.2226347 * c3 + 0.240390; 
      } 
      double x3 = x*x*x, x2 = x*x; 
      if(T >= 1667 && T <= 2222) { 
         y = -1.1063814 * x3 - 1.34811020 * x2 + 2.18555832 * x - 0.20219683; 
      } else if(T > 2222 && T <= 4000) { 
         y = -0.9549476 * x3 - 1.37418593 * x2 + 2.09137015 * x - 0.16748867; 
      } else if(T > 4000 && T <= 25000) { 
         y = 3.0817580 * x3 - 5.87338670 * x2 + 3.75112997 * x - 0.37001483; 
      }
      xyz[0] = x;
      xyz[1] = y;
      xyz[2] = 1 - x - y;
   }
}

double luma(double values[]) {
    double total = values[0] * 0.3086;
    total += values[1] * 0.6094;
    total += values[2] * 0.082;
    return total;
}

void tonemap(double rgb[], double avgloglum, double key, double burn) {
   double rhk = key / avgloglum;
   double rhb = (burn!=0) ? pow(2.0, burn) : 1.0;
   double Lp = luma(rgb)*rhk;
   scale(rgb, Lp * (1.0+Lp*rhb) / (1.0+Lp));
}

double log2(double x, double base) {
   return log(max(1e-37,x)) / log(max(1e-37,base));
}

double blackbody(double T, double adapt, double burn, double chr[]) {
    bbspline(T, chr);
    xyztochroma(chr);
    divide(chr, (C_PI * luma(chr)));
    double rgb[3] = {0,0,0};
    copy(chr, rgb);
    double k = max(1e-3,adapt);
    scale(rgb, 580.0 * log2(1.0 + planckpower(k*T), 1.0 + planckpower(k*5800)));
    tonemap(rgb,PLANCK_LPA_D10,0.18,burn);
    divide(rgb, PLANCK_TMNORM);
    return luma(rgb) / luma(chr);
}
