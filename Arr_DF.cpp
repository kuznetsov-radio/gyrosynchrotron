#include <math.h>
#include <float.h>
#include <stdlib.h>
#include "IDLinterface.h"
#include "Arr_DF.h"
#include "Plasma.h"

void Arr_DF :: Fp(double p, double p_z, double p_n, double *f, double *df_dp, double *df_dalpha)
{
 double G=sqrt(1.0+sqr(p/mc));
 double E=mc2*(G-1.0);

 if (iso_on)
 {
  double fE, dfE_dE;

  if (LQ_on) LQInterpolate(log_on ? log(E) : E, NE, E_arr, f_avg, &fE, &dfE_dE);
  else S1->Interpolate(log_on ? log(E) : E, &fE, &dfE_dE);

  if (log_on)
  {
   fE=exp(fE);
   dfE_dE*=(fE/E);
  }

  *f=fE/(p*me*G); 
  *df_dp=(dfE_dE-fE*G*me/sqr(p)*(1.0+sqr(p/G/mc)))/sqr(me*G);
  *df_dalpha=0.0;
 }
 else
 {
  double mu=(p>0.0) ? p_z/p : 0.0;
  if (mu>1.0) mu=1.0;
  if (mu<(-1.0)) mu=-1.0;
  double sa=(p>0.0) ? p_n/p : 1.0;

  double fE, dfE_dE, dfE_dmu;
  if (!LQ_on) S2->Interpolate(log_on ? log(E) : E, mu, &fE, &dfE_dE, &dfE_dmu, 0);
  else LQInterpolate2D(log_on ? log(E) : E, mu, NE, Nmu, E_arr, mu_arr, f_arr, &fE, &dfE_dE, &dfE_dmu, 0);

  if (log_on)
  {
   fE=exp(fE);
   dfE_dE*=(fE/E);
   dfE_dmu*=fE;
  }

  *f=fE/(p*me*G); 
  *df_dp=(dfE_dE-fE*G*me/sqr(p)*(1.0+sqr(p/G/mc)))/sqr(me*G);
  *df_dalpha=-dfE_dmu*sa/(p*me*G);
 }
}

void Arr_DF :: FE(double E, double mu, double *f, double *df_dE, double *df_dmu, double *g1, double *g2)
{
 if (iso_on)
 {
  double fE, dfE_dE;

  if (!f) *g1=0.0; //only g1 is needed
  else //calculating all
  {
   if (LQ_on) LQInterpolate(log_on ? log(E) : E, NE, E_arr, f_avg, &fE, &dfE_dE);
   else S1->Interpolate(log_on ? log(E) : E, &fE, &dfE_dE);

   if (log_on)
   {
    *f=exp(fE);
    *df_dE=dfE_dE*(*f/E);
   }
   else
   {
    *f=fE;
    *df_dE=dfE_dE;
   }

   *df_dmu=*g1=*g2=0.0;
  }
 }
 else
 {
  double fE, dfE_dE, dfE_dmu, d2fE_dmu2;

  if (!f) //only g1 is needed
  {
   if (LQ_on) LQInterpolate2D(log_on ? log(E) : E, mu, NE, Nmu, E_arr, mu_arr, f_arr, log_on ? 0 : &fE, 0, &dfE_dmu, 0);
   else S2->Interpolate(log_on ? log(E) : E, mu, log_on ? 0 : &fE, 0, &dfE_dmu, 0); 
   *g1=log_on ? dfE_dmu : dfE_dmu/fE; 
  }
  else //calculating all
  {
   if (LQ_on) LQInterpolate2D(log_on ? log(E) : E, mu, NE, Nmu, E_arr, mu_arr, f_arr, &fE, &dfE_dE, &dfE_dmu, &d2fE_dmu2);
   else S2->Interpolate(log_on ? log(E) : E, mu, &fE, &dfE_dE, &dfE_dmu, &d2fE_dmu2);

   if (log_on)
   {
    *f=exp(fE);
    *df_dE=dfE_dE*(*f/E);
    *df_dmu=dfE_dmu*(*f);
    *g1=dfE_dmu;
    *g2=d2fE_dmu2+sqr(*g1);
   }
   else
   {
    *f=fE;
    *df_dE=dfE_dE;
    *df_dmu=dfE_dmu;
    *g1=dfE_dmu/fE; 
    *g2=d2fE_dmu2/fE;
   }
  }
 }
}

double Arr_DF :: IntegratedF()
{
 double *a=(double*)malloc(sizeof(double)*Nmu);
 double *Eprof=(double*)malloc(sizeof(double)*NE);

 for (int j=0; j<Nmu; j++) 
 {
  for (int i=0; i<NE; i++) Eprof[i]=f_arr[i][j];
  a[j]=log_on ? IntTabulatedLog(E_arr, Eprof, NE) : IntTabulated(E_arr, Eprof, NE); 
 }
 double res=2.0*M_PI*IntTabulated(mu_arr, a, Nmu);

 free(a);
 free(Eprof);

 return res;
}

Arr_DF :: Arr_DF(int *Lparms, double *E, double *mu, double *f, int *OK, int *empty)
{
 E_arr=mu_arr=f_avg=0;
 f_arr=0;
 S1=0;
 S2=0;

 *OK=1;

 if (Lparms[i_NE]<3 || Lparms[i_Nmu]<3) *empty=1;
 else
 {
  NE=Lparms[i_NE];
  Nmu=Lparms[i_Nmu];

  log_on=Lparms[i_logkey]==0;

  iso_on=PK_on=0;
  if (Lparms[i_PKkey]==1) iso_on=1;
  if (Lparms[i_PKkey]==2) iso_on=PK_on=1;

  LQ_on=Lparms[i_splinekey]!=0;

  E_arr=(double*)malloc(sizeof(double)*NE);
  for (int i=0; i<NE; i++)
  {
   double q=E[i]*eV*1e6;
   if (log_on && q<=0)
   {
	*OK=0;
	break;
   }
   else E_arr[i]=log_on ? log(q) : q;
  }

  if (*OK)
  {
   for (int i=1; i<NE; i++) if (E_arr[i]<=E_arr[i-1])
   {
	*OK=0;
	break;
   }

   if (*OK)
   {
    N_intervals=1;
	E_x[0]=log_on ? exp(E_arr[0]) : E_arr[0];
	E_x[1]=log_on ? exp(E_arr[NE-1]) : E_arr[NE-1];
	logscale[0]=log_on;

    mu_arr=(double*)malloc(sizeof(double)*Nmu);
    for (int j=0; j<Nmu; j++) mu_arr[j]=mu[j];

	for (int j=1; j<Nmu; j++) if (mu_arr[j]<=mu_arr[j-1])
	{
	 *OK=0;
	 break;
	}

	if (*OK)
	{
	 if (iso_on) f_avg=(double*)malloc(sizeof(double)*NE);
	 double *mu_prof=(double*)malloc(sizeof(double)*Nmu);

	 f_arr=(double**)malloc(sizeof(double*)*NE);
	 for (int i=0; i<NE; i++) f_arr[i]=(double*)malloc(sizeof(double)*Nmu);

	 for (int i=0; i<NE; i++) 
	 {
	  for (int j=0; j<Nmu; j++)
	  {
	   double q=f[i+j*NE]/eV/1e6;
	   if (log_on && q<=0)
	   {
	    *OK=0; 
	    break;
	   }
	   else 
	   {
		f_arr[i][j]=log_on ? log(q) : q;
		mu_prof[j]=q;
	   }
	  }
	  if (!(*OK)) break;
	  if (iso_on) 
	  {
	   f_avg[i]=IntTabulated(mu_arr, mu_prof, Nmu)/(mu_arr[Nmu-1]-mu_arr[0]);
	   if (log_on) f_avg[i]=log(f_avg[i]);
	  }
	 }

	 free(mu_prof);

	 if (*OK)
	 {
	  nb=IntegratedF();
	  
	  if (!finite(nb) || (PosDef ? nb<0.0 : 0)) *OK=0;
	  *empty=(nb==0);

	  if (*OK && !(*empty))
	  {
	   if (!LQ_on) 
		if (iso_on) S1=new Spline(NE, E_arr, f_avg);
	    else S2=new Spline2D(NE, Nmu, E_arr, mu_arr, f_arr);

	   EPS_mu0=1e-3;

	   if (!iso_on)
	   {
	    double gmax=0.0;
        double g1;
        for (int i=0; i<NE; i++) for (int j=0; j<Nmu; j++)
        {
         FE(log_on ? exp(E_arr[i]) : E_arr[i], mu_arr[j], 0, 0, 0, &g1, 0);
         gmax=max(gmax, fabs(g1)); 
        }
        if (gmax!=0.0) EPS_mu0=min(EPS_mu0, 1.0/gmax/30); 
	   }
	  }
	 }
	}
   }
  }
 }
}

Arr_DF :: ~Arr_DF()
{
 if (E_arr) free(E_arr);
 if (mu_arr) free(mu_arr);
 if (f_avg) free(f_avg);
 if (f_arr)
 {
  for (int i=0; i<NE; i++) free(f_arr[i]);
  free(f_arr);
 }
 if (S1) delete S1;
 if (S2) delete S2;
}