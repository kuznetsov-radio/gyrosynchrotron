#include <math.h>
#include <float.h>
#include <malloc.h>
#include <time.h>
#include "ExtMath.h"
#include "IDLInterface.h"
#include "DF.h"
#include "GS.h"
#include "Astrophys.h"
#include "Messages.h"

#ifdef LINUX
#define CLK_TCK CLOCKS_PER_SEC
#endif

EmWave :: EmWave(double _nu, double _theta, int _sigma, double _nu_p, double _nu_B)
{
 nu=_nu;
 theta=_theta;
 nu_p=_nu_p;
 nu_B=_nu_B;
 sigma=_sigma;

 double nu_c=(sigma==-1) ? nu_B/2+sqrt(sqr(nu_p)+sqr(nu_B)/4) : nu_p; //cutoff frequency
 
 if (nu<=nu_c) Valid=0;
 else
 {
  Valid=1;

  ct=cos(theta);
  st=sin(theta);
  if (fabs(ct)<cst_min)
  {
   ct=cst_min*sign(ct);
   st=sqrt(1.0-sqr(ct))*sign(st);
  }
  if (fabs(st)<cst_min)
  {
   st=cst_min*sign(st);
   ct=sqrt(1.0-sqr(st))*sign(ct);
  }

  y=nu/nu_B;

  double u=sqr(nu_B/nu);
  double v=sqr(nu_p/nu);
  double Delta=sqrt(sqr(u*sqr(st))+4.0*u*sqr((1.0-v)*ct));
  N=sqrt(1.0-2.0*v*(1.0-v)/(2.0*(1.0-v)-u*sqr(st)+double(sigma)*Delta)); //refraction index
  N_z=N*ct; //longitudinal component of the refraction index
  T=2.0*sqrt(u)*(1.0-v)*ct/(u*sqr(st)-double(sigma)*Delta); //axial polarization coefficient;
  L=(v*sqrt(u)*st+T*u*v*st*ct)/(1.0-u-v+u*v*sqr(ct)); //longitudinal polarization coefficient
  Zfactor=u ? (2.0*sqr(1.0-v)+u*sqr(st)-double(sigma)*sqr(u*sqr(st))/Delta)/
	          sqr(2.0*(1.0-v)-u*sqr(st)+double(sigma)*Delta) : 0.5; 
              //Zheleznyakov's correction to free-free

  Valid=finite(N); 
 }
}

//--------------------------------------------------------------------------

typedef struct
{
 double Q, df_dp, df_dalpha;
} CacheStruct;

class GSIntegrand : public IntegrableFunction
{
 public:
 EmWave *w;
 DF *df;
 int ExactBessel;
 
 int s, mode;
 double x;
 
 CacheStruct *Cache;
 int rom_count, rom_n, CacheSize;

 double F(double p_z);
};

double GSIntegrand :: F(double p_z)
{
 double G=p_z*w->N_z/mc+x;
 double p2=sqr(mc)*(sqr(G)-1.0);
 double p=(p2>0.0) ? sqrt(p2) : 0.0;
 double p_n2=p2-sqr(p_z);   
 double p_n=(p_n2>0.0) ? sqrt(p_n2) : 0.0;
 double ca=(p2>0.0) ? p_z/p : 0.0;
 double sa=(p2>0.0) ? p_n/p : 1.0;
 double beta=p/G/mc;
 double beta_z=p_z/G/mc;
 double beta_n=p_n/G/mc;

 double Q, f, df_dp, df_dalpha;

 int i=rom_count++;

 if (rom_count>rom_n)
 {
  double lambda=w->y*p_n/mc*w->N*w->st;
  if (lambda!=0.0)
  {
   double Js, Js1;
   if (ExactBessel) FindBesselJ(lambda, s, &Js, &Js1);       
   else FindBesselJ_WH(lambda, s, &Js, &Js1);
   Q=sqr((w->T*(w->ct-w->N*beta_z)+w->L*w->st)/(w->N*beta_n*w->st)*Js+Js1);
  }
  else if (s==1) Q=sqr(((w->T*(w->ct-w->N*beta_z)+w->L*w->st)*G*w->y+1.0)/2);
  else Q=0.0;
  df->Fp(p, p_z, p_n, &f, &df_dp, &df_dalpha);

  if (i>=CacheSize)
  {
   CacheSize+=64;
   Cache=(CacheStruct*)realloc(Cache, CacheSize*sizeof(CacheStruct));
  }

  Cache[i].Q=Q;
  Cache[i].df_dp=df_dp;
  Cache[i].df_dalpha=df_dalpha;
  rom_n=rom_count;
 }
 else
 {
  Q=Cache[i].Q;
  df_dp=Cache[i].df_dp;
  df_dalpha=Cache[i].df_dalpha;
 }

 Q*=(mode ? G*mc*sa*(p_n*df_dp+(ca-w->N_z*beta)*df_dalpha) : f*sqr(p_n));

 return Q;
}

void GS_jk(EmWave *w, DF *df, int ExactBessel, double *j, double *k)
{
 if (!w->Valid)
 {
  *j=0.0;
  *k=1e100;
 }
 else
 {
  *j=*k=0.0;

  if (w->nu_B>0.0)
  {
   GSIntegrand gsi;
   gsi.ExactBessel=ExactBessel;
   gsi.w=w;
   gsi.df=df;

   gsi.CacheSize=65;
   gsi.Cache=(CacheStruct*)malloc(gsi.CacheSize*sizeof(CacheStruct));

   int timeout=0;
   int s_out=0;
   clock_t t0=clock();

   for (int r=0; r<df->N_intervals && !timeout && !s_out; r++)
   {
    double j_loc, k_loc;
    j_loc=k_loc=0.0;

    //calculating the distribution function parameters:
	double E_min=df->E_x[r];
	double E_max=df->E_x[r+1];
	if (r>0) E_min*=(1.0+1e-10);
	if (r<(df->N_intervals-1)) E_max*=(1.0-1e-10);
    double G_min=E_min/mc2+1.0; //minimal Lorentz factor
    double G_max=E_max/mc2+1.0; //maximal Lorentz factor
    double p_min=mc*sqrt(sqr(G_min)-1.0); //minimal impulse
    double p_max=mc*sqrt(sqr(G_max)-1.0); //maximal impulse

    int p_out=0;
    int j_done=0;
    int k_done=0;
    int vfinite=1;
      
    gsi.s=S_MIN;
    do
    {
     gsi.x=w->nu_B/w->nu*gsi.s;
     double R2=sqr(w->N_z)+(gsi.x-1.0)*(gsi.x+1.0);
     if (R2>0.0)
     {
      double R=sqrt(R2);
	  double p_z1=mc*(w->N_z*gsi.x-R)/(1.0-sqr(w->N_z));       
	  double p_z2=mc*(w->N_z*gsi.x+R)/(1.0-sqr(w->N_z));       

      int M=1;
      if (p_z1>p_max || p_z2<(-p_max) || (p_z1>(-p_min) && p_z2<p_min)) M=0;
      if (p_z1<(-p_max) && p_z2>p_max)
      {
       M=0;
       if (gsi.x>1) p_out=1;
      }

      if (M)
      {
       double pzx=mc/w->N_z*(G_min-gsi.x);
       if (fabs(pzx)<p_min)
       {
        if (w->N_z>0) p_z1=pzx;
        else p_z2=pzx;
       }
       pzx=mc/w->N_z*(G_max-gsi.x);
       if (fabs(pzx)<p_max)
       {
        if (w->N_z>0) p_z2=pzx;
        else p_z1=pzx;
       }

       int err;
       double q;
       gsi.rom_n=0;
       if (!j_done)
       {
        gsi.mode=0;
        gsi.rom_count=0;
        q=qromb(&gsi, p_z1, p_z2, ERR_i, &err);
        j_loc+=q;
        if (j_loc!=0.0) if (fabs(q/j_loc)<ERR_s) j_done=1;
	    if (!finite(q)) vfinite=0;
       }

       if (!k_done)
       {
        gsi.mode=1;
        gsi.rom_count=0;
        q=qromb(&gsi, p_z1, p_z2, ERR_i, &err); 
        k_loc+=q;
        if (k_loc!=0.0) if (fabs(q/k_loc)<ERR_s) k_done=1;
	    if (!finite(q)) vfinite=0;
       }
      }
     }

     if (gsi.s>=S_MAX) s_out=1;
     if (((clock()-t0)/CLK_TCK)>T_MAX) timeout=1;

     gsi.s++; 
    } while(!p_out && !(j_done && k_done) && !s_out && !timeout && vfinite);

    *j+=j_loc;
    *k+=k_loc;
   }

   *j*=(4.0*sqr(M_PI*e)/c*w->N*w->nu/(1.0+sqr(w->T)));
   *k*=(-4.0*sqr(M_PI*e)/w->N/w->nu/(1.0+sqr(w->T)));
 
   free(gsi.Cache);
  }
 }
}

//--------------------------------------------------------------------

typedef struct
{
 double Q, R;
} CacheStructApprox;

class GSIntegrandApprox : public IntegrableFunction
{
 public:
 EmWave *w;
 DF *df;
 
 int Q_on;
 int Qcorr, mode;
 double mu_list[2], lnQ1_list[2];
 int mflag;
 
 CacheStructApprox *Cache;
 int rom_count, rom_n, CacheSize;

 double E_loc, beta_loc, G_loc, p_loc;
 double H1(double mu);
 void Find_mu0(double *mu0, double *lnQ2);

 double F(double E);
};

double GSIntegrandApprox :: H1(double mu)
{
 double g1;
 df->FE(E_loc, mu, 0, 0, 0, &g1, 0);

 double nbct=w->N*beta_loc*w->ct;
 double nbmct1=1.0-nbct*mu;
 double sa2=1.0-sqr(mu);
 double sa=sqrt(sa2);
 double x=w->N*beta_loc*w->st*sa/nbmct1;
 double s1mx2=sqrt(1.0-sqr(x));
 double lnZ=s1mx2+log(x/(1.0+s1mx2));

 double lnQ1=0.0;
 if (Qcorr)
 {
  double s1mx2_3=s1mx2*s1mx2*s1mx2;
  double s=G_loc*w->y*nbmct1;
  double a6=s1mx2_3+0.503297/s;
  double b16=s1mx2_3+1.193000/s;
  double b2=(1.0-0.2*pow(s, -2.0/3));
  double ab=pow(a6*b16, 1.0/6)*b2;
  double xi=3.0*sqr(x)*s1mx2*(w->N_z*beta_loc-mu)/sa2;
  double eta=w->N_z*beta_loc/s;
  double lambda=G_loc*w->y/(6.0*s);
  double a_1a=lambda*(0.503297*eta-xi)/a6;
  double b_1b=lambda*(1.193000*eta-xi)/b16+4.0*lambda*beta_loc*w->N_z*(b2-1.0)/b2;
  lnQ1=2.0*(ab*(a_1a+b_1b)*nbmct1-ab*w->N_z*beta_loc-w->T*w->N*beta_loc)/
       (w->T*(w->ct-w->N*beta_loc*mu)+w->L*w->st+ab*nbmct1)-2.0*a_1a+w->N_z*beta_loc/nbmct1;

  mu_list[mflag]=mu;
  lnQ1_list[mflag]=lnQ1;
  mflag^=1;
 }

 return g1+2.0*G_loc*w->y*((nbct-mu)/sa2*s1mx2-nbct*lnZ)+lnQ1;
}

class MuSolveFunction : public IntegrableFunction
{
 public:
 GSIntegrandApprox *gsi;
 double F(double mu);
};

double MuSolveFunction :: F(double mu)
{ 
 double h=gsi->H1(mu); 
 return h;
}

void GSIntegrandApprox :: Find_mu0(double *mu0, double *lnQ2)
{
 MuSolveFunction msf;
 msf.gsi=this;

 *lnQ2=0.0;

 Qcorr=0;
 *mu0=BrentRoot(&msf, -1.0+1e-5, 1.0-1e-5, df->EPS_mu0);

 if (finite(*mu0) && Q_on)
 {
  mflag=0;
  int Qfound=0;

  Qcorr=1;
  double mu1=SecantRoot(&msf, *mu0, *mu0-1e-4*sign(*mu0), df->EPS_mu0);
 
  if (finite(mu1) && fabs(mu1)<1.0) 
  {
   *mu0=mu1;
   Qfound=1;   
  }
  else
  {
   mu1=BrentRoot(&msf, -1.0+1e-5, *mu0, df->EPS_mu0);
   if (finite(mu1))
   {
	*mu0=mu1;
	Qfound=1;
   }
   else
   {
	mu1=BrentRoot(&msf, *mu0, 1.0-1e-5, df->EPS_mu0);
	if (finite(mu1))
	{
	 *mu0=mu1;
	 Qfound=1;
	}
   }
  }

  if (Qfound) *lnQ2=(lnQ1_list[1]-lnQ1_list[0])/(mu_list[1]-mu_list[0]);
 }
}

double GSIntegrandApprox :: F(double E)
{
 if (E==0.0) return 0.0;
 else
 {
  double f, Q, R;

  int i=rom_count++;

  if (rom_count>rom_n)
  {
   E_loc=E;
   G_loc=E/mc2+1.0;
   beta_loc=sqrt(sqr(G_loc)-1.0)/G_loc;
   p_loc=beta_loc*G_loc*mc;

   double mu0, lnQ2;
   if (df->PK_on)
   {
    mu0=w->N*beta_loc*w->ct;
	lnQ2=0;
   }
   else Find_mu0(&mu0, &lnQ2);

   if (finite(mu0))
   {
	double df_dE, df_dmu, g1, g2;
    df->FE(E, mu0, &f, &df_dE, &df_dmu, &g1, &g2);

	double nbct=w->N*beta_loc*w->ct;
	double nbctm=nbct-mu0;
    double nbmct1=1.0-nbct*mu0;
    double sa2=1.0-sqr(mu0);
    double sa=sqrt(sa2);
    double x=w->N*beta_loc*w->st*sa/nbmct1;
    double s1mx2=sqrt(1.0-sqr(x));
	double s1mx2_3=s1mx2*s1mx2*s1mx2;
	double s=G_loc*w->y*nbmct1;
	double a=pow(s1mx2_3+0.503297/s, 1.0/6);
	double b=pow(s1mx2_3+1.193000/s, 1.0/6)*(1.0-0.2*pow(s, -2.0/3));

	Q=sqr(w->T*(w->ct-w->N*beta_loc*mu0)+w->L*w->st+a*b*nbmct1)/sqr(a)/nbmct1;

	double lnZ=s1mx2+log(x/(1.0+s1mx2));
	double Z=exp(2.0*s*lnZ);

	double H2=g2-sqr(g1)-2.0*G_loc*w->y*s1mx2/sa2*(1.0+sqr(w->N*beta_loc*w->st*nbctm/nbmct1)/nbmct1/(1.0-sqr(x))-
		      2.0*mu0*nbctm/sa2+nbct*nbctm/nbmct1)+lnQ2;
	double LpFactor=sqrt(-2.0*M_PI/H2);
	if (!finite(LpFactor)) LpFactor=0;

	Q*=(Z*LpFactor);

	R=df_dE-(1.0+sqr(beta_loc))/(c*p_loc*beta_loc)*f+nbctm/(c*p_loc*beta_loc)*df_dmu;
   } 
   else f=Q=R=0.0; 

   if (i>=CacheSize)
   {
    CacheSize+=64;
    Cache=(CacheStructApprox*)realloc(Cache, CacheSize*sizeof(CacheStructApprox));
   }

   Cache[i].Q=Q;
   Cache[i].R=R;
   rom_n=rom_count;
  }
  else
  {
   Q=Cache[i].Q;
   R=Cache[i].R;
  }

  return mode ? Q*R : Q*f;
 }
}

void GS_jk_approx(EmWave *w, DF *df, int Npoints, int Q_on, double *j, double *k)
{
 if (!w->Valid)
 {
  *j=0.0;
  *k=1e100;
 }
 else
 {
  *j=*k=0.0;

  if (w->nu_B>0.0)
  {
   GSIntegrandApprox gsi;
   gsi.w=w;
   gsi.df=df;
   gsi.Q_on=Q_on;

   gsi.CacheSize=65;
   gsi.Cache=(CacheStructApprox*)malloc(gsi.CacheSize*sizeof(CacheStructApprox));

   for (int r=0; r<df->N_intervals; r++)
   {
    double jk_loc[2];
    int err;
    gsi.rom_n=0;

    for (gsi.mode=0; gsi.mode<=1; gsi.mode++)
    {
     gsi.rom_count=0; 
     jk_loc[gsi.mode]=df->logscale[r] ? 
	 	              ((Npoints>=1) ? trapzdLog(&gsi, df->E_x[r], df->E_x[r+1], Npoints) : 
					                  qrombLog( &gsi, df->E_x[r], df->E_x[r+1], ERR_i, &err)) : 
		              ((Npoints>=1) ? trapzd(   &gsi, df->E_x[r], df->E_x[r+1], Npoints) : 
						              qromb(    &gsi, df->E_x[r], df->E_x[r+1], ERR_i, &err)); 
    }
 
    *j+=jk_loc[0];
    *k+=jk_loc[1];
   }

   *j*=(2.0*M_PI*sqr(e)*w->nu/c/w->N/(1.0+sqr(w->T))/sqr(w->st));
   *k*=(-2.0*M_PI*sqr(e)*c/sqr(w->N)/w->N/w->nu/(1.0+sqr(w->T))/sqr(w->st));

   free(gsi.Cache);
  }
 }
}

//--------------------------------------------------------------

void GS_jk_mDF(EmWave *w, DF **df, int ExactBessel, double *j, double *k)
{
 *j=*k=0;

 DF **df_loc=df;

 while (*df_loc)
 {
  double j_loc, k_loc;

  GS_jk(w, *df_loc, ExactBessel, &j_loc, &k_loc);

  *j+=j_loc;
  *k+=k_loc;

  df_loc++;
 }
}

void GS_jk_approx_mDF(EmWave *w, DF **df, int Npoints, int Q_on, double *j, double *k)
{
 *j=*k=0;

 DF **df_loc=df;

 while (*df_loc)
 {
  double j_loc, k_loc;

  GS_jk_approx(w, *df_loc, Npoints, Q_on, &j_loc, &k_loc);
                  
  *j+=j_loc;
  *k+=k_loc;

  df_loc++;
 }
}

//----------------------------------------------------

double Coulomb(double T0, double nu)
{
 return (T0<2e5) ? 18.2+1.5*log(T0)-log(nu) : 24.573+log(T0/nu);
}

inline double Zmean(double T0)
{
 return (T0>3.5e4) ? 1.146 : 1.0;
}

double Hminus(EmWave *w, double ne, double nH, double T0)
{
 /* Eq. (20) from Stallcop (ApJ, 187, 179, 1974).
    We consider this effect only if T<160e3 K which corresponds to kT<1; 
    in the original paper, it was analyzed only for T<50e3 K. */
 double Hm=0;
 if (T0<160e3)
 {
  double kT=sqrt(kB*T0/ieH);
  double xi=4.862*kT*(1.0-0.2096*kT+0.0170*kT*kT-0.00968*kT*kT*kT); 
  Hm=nH ? 2.15e-29*ne*nH*kB*T0*sqr(ieH/hPl/w->nu)*exp(-xi)/kT/w->N : 0;
 }
 return Hm;
}

void FF_jk_Maxwell(EmWave *w, double ne, double nH, double T0, double *j, double *k)
{
 if (!w->Valid)
 {
  *j=0.0;
  *k=1e100;
 }
 else if (ne==0.0)
 {
  *k=0.0;
  *j=0.0;
 }
 else
 {
  double lnL=Coulomb(T0, w->nu);
  double kFF=9.786e-3*sqr(ne)*lnL/(w->N*sqr(w->nu)*T0*sqrt(T0));
  double kH=Hminus(w, ne, nH, T0);

  *k=(kFF+kH)*2.0*w->Zfactor*Zmean(T0); 
  *j=(*k)*sqr(w->N*w->nu/c)*kB*T0; 
 }
}

void FF_jk_kappa(EmWave *w, double ne, double nH, double T0, double kappa, double *j, double *k)
{
 if (!w->Valid)
 {
  *j=0.0;
  *k=1e100;
 }
 else if (ne==0.0)
 {
  *j=0.0;
  *k=0.0;
 }
 else
 {
  double lnL=Coulomb(T0, w->nu);

  double Ak=Gamma(kappa+1.0)/Gamma(kappa-0.5)*pow(kappa-1.5, -1.5);
  double kFF=Ak*8.0*sqr(e*e*e)*sqr(ne)*lnL/
	         (3.0*sqrt(2.0*M_PI)*w->N*c*sqr(w->nu)*m*kB*T0*sqrt(m*kB*T0))*
			 (1.0-0.575*pow(6.0/kappa, 1.1)/lnL);
  double jFF=Ak*(kappa-1.5)/kappa*8.0*sqr(e*e*e)*w->N*sqr(ne)*lnL/
	         (3.0*sqrt(2.0*M_PI)*mc2*sqrt(mc2)*sqrt(kB*T0))*
			 (1.0-0.525*pow(4.0/kappa, 1.25)/lnL);
  kFF*=(2.0*w->Zfactor*Zmean(T0));
  jFF*=(2.0*w->Zfactor*Zmean(T0));

  double kH=Hminus(w, ne, nH, T0);
  kH*=(2.0*w->Zfactor*Zmean(T0));
  double jH=kH*sqr(w->N*w->nu/c)*kB*T0;

  *k=kFF+kH; 
  *j=jFF+jH;
 }
}

void FF_jk(EmWave *w, double ne, double nH, double T0, double kappa, double *j, double *k)
{
 if (!finite(kappa)) FF_jk_Maxwell(w, ne, nH, T0, j, k);
 else FF_jk_kappa(w, ne, nH, T0, kappa, j, k);
}
