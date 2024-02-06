#include <stdlib.h>
#include <math.h>
#include "IDLinterface.h"
#include "Messages.h"
#include "Plasma.h"
#include "ExtMath.h"
#include "DF.h"
#include "Std_DF.h"
#include "Arr_DF.h"
#include "GS.h"
#include "FF.h"

int FindLocalJK(double *nu, int *Lparms, double *Rparms, double *Parms,  
	            double *E_arr, double *mu_arr, double *f_arr, 
	            double *jX, double *jO, double *kX, double *kO, double *ne_total) 
{
 int res=0;

 int EM_flag=int(Parms[i_EMflag]);
 int DFtypeGlobalFlag=Lparms[i_arrKeyG];
 int DFtypeLocalFlag=int(Parms[i_arrKeyL]);
 double kappa=dNaN;

 int Ndf=0;
 DF *df[10];
 
 //initializing the analytical distribution function:
 if (!(DFtypeGlobalFlag & 2)) if (!(DFtypeLocalFlag & 2))
 {
  int k=0;
  int Done=0;
  int OK, empty, kap_on;

  while (!Done)
  {
   df[Ndf]=new Std_DF(Parms, k, &OK, &empty, &kap_on, &Done);
   if (!OK) res=1;

   if (OK && !empty) 
   {
    Ndf++;
    if (kap_on) kappa=Parms[i_epskappa];
   }
   else delete df[Ndf];
   k++;
  }
 }

 //initializing the array distribution function
 if (!(DFtypeGlobalFlag & 1)) if (!(DFtypeLocalFlag & 1))
 {
  int OK, empty;

  df[Ndf]=new Arr_DF(Lparms, E_arr, mu_arr, f_arr, &OK, &empty);
  
  if (!OK) res=2;

  if (OK && !empty) Ndf++;
  else delete df[Ndf];
 }
 
 df[Ndf]=0;
 
 if (!res)
 {
  double nb=0.0; //additional energetic electron density
  for (int i=0; i<Ndf; i++) nb+=df[i]->nb;

  *ne_total=Parms[i_n0]+nb;
  
  double nu_p=e*sqrt(*ne_total/me/M_PI); 
  double nu_B=e*Parms[i_B]/me/c/(2.0*M_PI); 
  double theta=Parms[i_theta]*M_PI/180;

  int Nnu=Lparms[i_Nnu];
  for (int i=0; i<Nnu; i++) jX[i]=kX[i]=jO[i]=kO[i]=0;

  if (!(EM_flag & 1)) if (Ndf) //GS is on and nonthermal electrons are present
  {
   double nu_cr=nu_B*Rparms[i_nuCr];
   double nu_cr_WH=nu_B*Rparms[i_nuWH];

   int Npoints=Lparms[i_Nnodes];
   if (Npoints>=0) Npoints=max(Npoints, 16);
   if (Npoints<0) Npoints=0;

   int Q_on=Lparms[i_QoptKey]==0;
   int m_on=Lparms[i_matchKey]==0;

   Find_jk_GS(df, nu, Nnu, -1, theta, nu_p, nu_B, nu_cr, nu_cr_WH, Npoints, Q_on, m_on, jX, kX);
   Find_jk_GS(df, nu, Nnu,  1, theta, nu_p, nu_B, nu_cr, nu_cr_WH, Npoints, Q_on, m_on, jO, kO);
  }

  if (!(EM_flag & 2)) //e-ions is on
  {
   double j_loc, k_loc;

   for (int i=0; i<Nnu; i++)
   {
	Find_jk_FFei(Parms[i_n0], Parms[i_T0], nu_p, nu_B, theta, kappa, int(Parms[i_abcode]), -1, nu[i], &j_loc, &k_loc);
	jX[i]+=j_loc;
	kX[i]+=k_loc;
	Find_jk_FFei(Parms[i_n0], Parms[i_T0], nu_p, nu_B, theta, kappa, int(Parms[i_abcode]),  1, nu[i], &j_loc, &k_loc);
	jO[i]+=j_loc;
	kO[i]+=k_loc;
   }
  }

  if (!(EM_flag & 4)) //e-neutrals is on
  {
   double j_loc, k_loc;

   for (int i=0; i<Nnu; i++)
   {
	Find_jk_FFen(Parms[i_n0], Parms[i_nH], Parms[i_nHe], Parms[i_T0], nu_p, nu_B, theta, -1, nu[i], &j_loc, &k_loc);
	jX[i]+=j_loc;
	kX[i]+=k_loc;
	Find_jk_FFen(Parms[i_n0], Parms[i_nH], Parms[i_nHe], Parms[i_T0], nu_p, nu_B, theta,  1, nu[i], &j_loc, &k_loc);
	jO[i]+=j_loc;
	kO[i]+=k_loc;
   }
  }
 }

 for (int i=0; i<Ndf; i++) delete(df[i]);

 return res;
}

void RadiationTransfer(double nu, int Nz, double *dz, double *ne, double *B, double *theta, 
	                   double *jX, double *jO, double *kX, double *kO, 
					   double *Lw, double *Rw, double *Ls, double *Rs, double *Le, double *Re)
{
 for (int i=0; i<Nz; i++)
 {
  double tau=-kO[i]*dz[i];
  double eO=(tau<700) ? exp(tau) : 0.0; 
  double dIO=(kO[i]==0.0 || tau>700) ? 0.0 : jO[i]/kO[i]*((1.0-eO) ? 1.0-eO : -tau);
  tau=-kX[i]*dz[i];
  double eX=(tau<700) ? exp(tau) : 0.0;
  double dIX=(kX[i]==0.0 || tau>700) ? 0.0 : jX[i]/kX[i]*((1.0-eX) ? 1.0-eX : -tau);

  if (i>0) if (((theta[i]>(M_PI/2)) ^ (theta[i-1]>(M_PI/2))) && ne[i]>0 && ne[i-1]>0)
  {
   double a=*Lw;
   *Lw=*Rw;
   *Rw=a;
   
   double B_avg=(B[i]+B[i-1])/2;
   double ne_avg=(ne[i]+ne[i-1])/2;
   double da_dz=fabs(theta[i]-theta[i-1])/(dz[i]+dz[i-1])*2;

   double QT=e*e*e*e*e/(32*M_PI*M_PI*me*me*me*me*c*c*c*c)*ne_avg*sqr(B_avg)*B_avg/sqr(sqr(nu))/da_dz;
   QT=exp(-QT);
   a=*Le*QT+*Re*(1.0-QT);
   *Re=*Re*QT+*Le*(1.0-QT);
   *Le=a;
  }

  if (theta[i]>(M_PI/2))
  {
   *Lw=dIX+*Lw*eX;
   *Ls=dIX+*Ls*eX;
   *Le=dIX+*Le*eX;
   *Rw=dIO+*Rw*eO;
   *Rs=dIO+*Rs*eO;
   *Re=dIO+*Re*eO;
  }
  else
  {
   *Lw=dIO+*Lw*eO;
   *Ls=dIO+*Ls*eO;
   *Le=dIO+*Le*eO;
   *Rw=dIX+*Rw*eX;
   *Rs=dIX+*Rs*eX;
   *Re=dIX+*Re*eX;
  }
 }
}

int MW_Transfer(int *Lparms, double *Rparms, double *Parms, double *E_arr, double *mu_arr, double *f_arr, double *RL)
{
 int res=0;

 int Nnu=Lparms[i_Nnu];
 double *nu=(double*)malloc(sizeof(double)*Nnu);
 if (Rparms[i_nu0]>0)
 {
  nu[0]=Rparms[i_nu0];
  double dnu=pow(10.0, Rparms[i_dnu]);
  for (int i=1; i<Nnu; i++) nu[i]=nu[i-1]*dnu;
 }
 else for (int i=0; i<Nnu; i++) nu[i]=RL[i*OutSize+iRL_nu]*1e9;

 int Nz=Lparms[i_Nz];
 double **jX, **jO, **kX, **kO, *dz, *ne_total, *B, *theta;
 jX=(double**)malloc(sizeof(double*)*Nz);
 jO=(double**)malloc(sizeof(double*)*Nz);
 kX=(double**)malloc(sizeof(double*)*Nz);
 kO=(double**)malloc(sizeof(double*)*Nz);
 dz=(double*)malloc(sizeof(double)*Nz);
 ne_total=(double*)malloc(sizeof(double)*Nz);
 B=(double*)malloc(sizeof(double)*Nz);
 theta=(double*)malloc(sizeof(double)*Nz);
 for (int i=0; i<Nz; i++)
 {
  jX[i]=(double*)malloc(sizeof(double)*Nnu);
  jO[i]=(double*)malloc(sizeof(double)*Nnu);
  kX[i]=(double*)malloc(sizeof(double)*Nnu);
  kO[i]=(double*)malloc(sizeof(double)*Nnu);
 }

 int err=0;
 for (int i=0; i<Nz && !err; i++)
 {
  dz[i]=Parms[i*InSize+i_dz];
  B[i]=Parms[i*InSize+i_B];
  theta[i]=Parms[i*InSize+i_theta]*M_PI/180;

  err=FindLocalJK(nu, Lparms, Rparms, Parms+i*InSize,  
	              E_arr, mu_arr, f_arr+i*Lparms[i_NE]*Lparms[i_Nmu], 
	              jX[i], jO[i], kX[i], kO[i], ne_total+i);
 }

 if (err) res=err;
 else
 {
  double Sang=Rparms[i_S]/(sqr(Rparms[i_dSun]*AU)*sfu);

  double Lw, Rw, Ls, Rs, Le, Re;

  double *jX_loc=(double*)malloc(sizeof(double)*Nz);
  double *jO_loc=(double*)malloc(sizeof(double)*Nz);
  double *kX_loc=(double*)malloc(sizeof(double)*Nz);
  double *kO_loc=(double*)malloc(sizeof(double)*Nz);

  for (int i=0; i<Nnu; i++) 
  {
   Lw=RL[i*OutSize+iRL_Lw]/Sang;
   Rw=RL[i*OutSize+iRL_Rw]/Sang;
   Ls=RL[i*OutSize+iRL_Ls]/Sang;
   Rs=RL[i*OutSize+iRL_Rs]/Sang;
   Le=RL[i*OutSize+iRL_Le]/Sang;
   Re=RL[i*OutSize+iRL_Re]/Sang;

   for (int j=0; j<Nz; j++)
   {
	jX_loc[j]=jX[j][i];
	jO_loc[j]=jO[j][i];
	kX_loc[j]=kX[j][i];
	kO_loc[j]=kO[j][i];
   }

   RadiationTransfer(nu[i], Nz, dz, ne_total, B, theta, jX_loc, jO_loc, kX_loc, kO_loc, &Lw, &Rw, &Ls, &Rs, &Le, &Re);

   RL[i*OutSize+iRL_nu]=nu[i]/1e9;
   RL[i*OutSize+iRL_Lw]=Lw*Sang;
   RL[i*OutSize+iRL_Rw]=Rw*Sang;
   RL[i*OutSize+iRL_Ls]=Ls*Sang;
   RL[i*OutSize+iRL_Rs]=Rs*Sang;
   RL[i*OutSize+iRL_Le]=Le*Sang;
   RL[i*OutSize+iRL_Re]=Re*Sang;
  }

  free(jX_loc);
  free(jO_loc);
  free(kX_loc);
  free(kO_loc);
 }

 for (int i=0; i<Nz; i++)
 {
  free(jX[i]);
  free(jO[i]);
  free(kX[i]);
  free(kO[i]);
 }
 free(jX);
 free(jO);
 free(kX);
 free(kO);
 free(dz);
 free(ne_total);
 free(B);
 free(theta);

 free(nu);

 return res;
}