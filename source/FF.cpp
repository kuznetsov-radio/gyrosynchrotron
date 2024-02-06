#include <math.h>
#include <float.h>
#include "ExtMath.h"
#include "Plasma.h"
#include "Coulomb.h"
#include "Zeta.h"

inline double CoulombOld(double T0, double nu)
{
 return (T0<2e5) ? 18.2+1.5*log(T0)-log(nu) : 24.573+log(T0/nu);
}

inline double ZmeanOld(double T0)
{
 return (T0>3.5e4) ? 1.146 : 1.0;
}

void FF_jk_Maxwell(EmWave *w, double ne, double T0, int ab, double *j, double *k)
{
 if (ab<0) //classical Dulk's formula
 {
  double lnC=CoulombOld(T0, w->nu);
  double kFF=9.786e-3*sqr(ne)*lnC/(w->N*sqr(w->nu)*T0*sqrt(T0));

  *k=kFF*w->Zfactor*ZmeanOld(T0);
  *j=(*k)*sqr(w->N*w->nu/c)*kB*T0; 
 }
 else //modern formulae with correct Coulomb logarithm and zeta-function
 {
  double lnC=lnC1(T0, w->nu);
  double zeta=Zeta_Solar(T0, w->nu, ab);

  double jff=8*e*e*e*e*e*e*w->N/(3.0*sqrt(2.0*M_PI)*sqrt(me*me*me)*c*c*c)*
             sqr(ne)*lnC/sqrt(kB*T0)*(1.0+zeta);
  double kff=8*e*e*e*e*e*e/(3.0*sqrt(2.0*M_PI)*w->N*c*sqr(w->nu)*sqrt(me*me*me))*
             sqr(ne)*lnC/(sqrt(kB*T0)*kB*T0)*(1.0+zeta);

  *j=jff*w->Zfactor;
  *k=kff*w->Zfactor;
 }
}

void FF_jk_kappa(EmWave *w, double ne, double T0, double kappa, int ab, double *j, double *k) //###
{
 double Ak=Gamma(kappa+1.0)/Gamma(kappa-0.5)*pow(kappa-1.5, -1.5);

 if (ab<0) //classical formulae from Fleishman & Kuznetsov (2014)
 {
  double lnC=CoulombOld(T0, w->nu);

  double kFF=Ak*8.0*sqr(e*e*e)*sqr(ne)*lnC/
	         (3.0*sqrt(2.0*M_PI)*w->N*c*sqr(w->nu)*me*kB*T0*sqrt(me*kB*T0))*
			 (1.0-0.575*pow(6.0/kappa, 1.1)/lnC);
  double jFF=Ak*(kappa-1.5)/kappa*8.0*sqr(e*e*e)*w->N*sqr(ne)*lnC/
	         (3.0*sqrt(2.0*M_PI)*mc2*sqrt(mc2)*sqrt(kB*T0))*
			 (1.0-0.525*pow(4.0/kappa, 1.25)/lnC);

  *k=kFF*w->Zfactor*ZmeanOld(T0);
  *j=jFF*w->Zfactor*ZmeanOld(T0);
 }
 else //more accurate formulae
 {
  double lnC=lnC1(T0, w->nu);
  double zeta=Zeta_Solar(T0, w->nu, ab);

  double kFF=Ak*8*e*e*e*e*e*e*sqr(ne)*lnC*(1.0+zeta)/
	         (3.0*sqrt(2.0*M_PI)*w->N*c*sqr(w->nu)*me*kB*T0*sqrt(me*kB*T0));
  double jFF=Ak*(kappa-1.5)/kappa*8*e*e*e*e*e*e*w->N*sqr(ne)*lnC*(1.0+zeta)/
	         (3.0*sqrt(2.0*M_PI)*mc2*sqrt(mc2)*sqrt(kB*T0));

  *k=kFF*w->Zfactor;
  *j=jFF*w->Zfactor;
 }
}

void Find_jk_FFei(double ne, double T0, double nu_p, double nu_B, double theta, double kappa, int abcode, 
	              int sigma, double nu, double *j, double *k) 
{
 EmWave w=EmWave(nu, theta, sigma, nu_p, nu_B, 0, 1);

 if (!w.Valid)
 {
  *j=0;
  *k=1e100;
 }
 else if (ne==0.0)
 {
  *j=0.0;
  *k=0.0;
 }
 else
 {
  int ab=0;
  if (abcode==1) ab=1;
  if (abcode==-1) ab=-1;

  if (finite(kappa)) FF_jk_kappa(&w, ne, T0, kappa, ab, j, k);
  else FF_jk_Maxwell(&w, ne, T0, ab, j, k);
 }
}

void Find_jk_FFen(double ne, double nH, double nHe, double T0, double nu_p, double nu_B, double theta, 
	              int sigma, double nu, double *j, double *k)
{
 EmWave w=EmWave(nu, theta, sigma, nu_p, nu_B, 0, 1);

 if (!w.Valid)
 {
  *j=0;
  *k=1e100;
 }
 else if (ne==0.0)
 {
  *j=0.0;
  *k=0.0;
 }
 else
 {
  double jH, kH;
  jH=kH=0;

  if (ne>0 && nH>0 && T0>2500 && T0<50000)
  {
   double kT=sqrt(kB*T0/ieH);
   double xi=4.862*kT*(1.0-0.2096*kT+0.0170*kT*kT-0.00968*kT*kT*kT);
   kH=1.2737207e-11*ne*nH*sqrt(T0)/(sqr(nu)*w.N)*exp(-xi);

   jH=kH*sqr(w.N*nu/c)*kB*T0;
  }

  double jHe, kHe;
  jHe=kHe=0;

  if (ne>0 && nHe>0 && T0>2500 && T0<25000)
  {
   double kT=sqrt(kB*T0/ieH);
   kHe=5.9375453e-13*ne*nHe*sqrt(T0)/(sqr(nu)*w.N)*(1.868+7.415*kT-22.56*kT*kT+15.59*kT*kT*kT);

   jHe=kHe*sqr(w.N*nu/c)*kB*T0;
  }

  *j=(jH+jHe)*w.Zfactor;
  *k=(kH+kHe)*w.Zfactor;
 }
}