#include <math.h>
#include "ExtMath.h"
#include "Astrophys.h"

double Saha(double n0, double T0)
{
 double x=0.0;
 if (T0>0.0)
 {
  double xi=pow(2.0*M_PI*m*kB*T0/sqr(hPl), 1.5)/n0*exp(-ieH/kB/T0); 
  x=xi ? 2.0/(sqrt(1.0+4.0/xi)+1.0) : 0.0;
 } 
 return x;
}