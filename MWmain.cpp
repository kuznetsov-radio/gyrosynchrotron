#include <malloc.h>
#include <math.h>
#include "ExtMath.h"
#include "IDLInterface.h"
#include "Astrophys.h"
#include "DF.h"
#include "GS.h"
#include "Messages.h"

void Find_jk_GS(DF **df, double *nu, int Nnu, int mode, double theta, double nu_p, double nu_B, double nu_cr, double nu_cr_WH, int Npoints, int Q_on, int m_on, double *j, double *k)
{
 for (int i=0; i<Nnu; i++)
 {
  EmWave w=EmWave(nu[i], theta, mode, nu_p, nu_B);
  if (nu[i]>nu_cr) GS_jk_approx_mDF(&w, df, Npoints, Q_on, j+i, k+i);
  else GS_jk_mDF(&w, df, nu[i]<nu_cr_WH, j+i, k+i);
 }

 if (m_on)
 {
  if (nu_cr!=nu_cr_WH && nu[0]<nu_cr && nu[Nnu-1]>nu_cr)
  {
   int i0;
   for (i0=1; i0<Nnu; i0++) if (nu[i0]>nu_cr) break;
 
   if (j[i0]!=0.0 && k[i0]!=0.0)
   {
    EmWave w=EmWave(nu[i0], theta, mode, nu_p, nu_B);
    double j0, k0;
	GS_jk_mDF(&w, df, nu[i0]<nu_cr_WH, &j0, &k0);
	double mj=j0/j[i0];
	double mk=k0/k[i0];
	for (int i=i0; i<Nnu; i++)
	{
	 j[i]*=mj;
	 k[i]*=mk; 
	}
   }
  }

  if (nu_cr_WH<=nu_cr && nu[0]<nu_cr_WH && nu[Nnu-1]>nu_cr_WH)
  {
   int i0;
   for (i0=1; i0<Nnu; i0++) if (nu[i0]>nu_cr_WH) break;
 
   if (j[i0]!=0.0 && k[i0]!=0.0)
   {
    EmWave w=EmWave(nu[i0], theta, mode, nu_p, nu_B);
    double j0, k0;
	GS_jk_mDF(&w, df, 1, &j0, &k0);
	double mj=j0/j[i0];
	double mk=k0/k[i0];
	for (int i=i0; i<Nnu; i++)
	{
	 j[i]*=mj;
	 k[i]*=mk; 
	}
   }
  }
 }

 if (nu_cr<0.0 && nu_cr_WH>0.0)
 {
  EmWave w=EmWave(nu_cr_WH, theta, mode, nu_p, nu_B);
  if (w.Valid)
  {
   double ja, ka;
   GS_jk_approx_mDF(&w, df, Npoints, Q_on, &ja, &ka);
   if (ja!=0.0 && ka!=0.0)
   {
    double je, ke;
    GS_jk_mDF(&w, df, 1, &je, &ke);
	double mj=je/ja;
	double mk=ke/ka;
	for (int i=0; i<Nnu; i++)
	{
	 j[i]*=mj;
	 k[i]*=mk;
	}
   }
  }
 }
}

void Add_jk_FF(double ne, double nH, double T0, double kappa, double *nu, int Nnu, int mode, double theta, double nu_p, double nu_B, double *j, double *k)
{
 if (T0>0) for (int i=0; i<Nnu; i++)
 {
  EmWave w=EmWave(nu[i], theta, mode, nu_p, nu_B);
  double jFF, kFF;
  FF_jk(&w, ne, nH, T0, kappa, &jFF, &kFF);
  j[i]+=jFF;
  k[i]+=kFF;
 }
}

int FindLocalJK(double *nu, int Nnu, double *ParmIn, int InSize, 
	            int NE, int Nmu, double *E_arr, double *mu_arr, IDL_DFarr *f3D,
	            double *jX, double *jO, double *kX, double *kO, double *ne_total)
{
 int res=0;

 double T0=ParmIn[2];
 
 double nH=D1(ParmIn, InSize, 29); //neutral hydrogen concentration
 double ne=D1(ParmIn, InSize, 30); //electron concentration (assumed to equal the proton concentration)
 if (nH+ne<=0) //if not defined directly, calculate using Saha equation
 {
  double n0=ParmIn[11];
  double alpha=Saha(n0, fabs(T0));
  ne=n0*alpha;
  nH=n0*(1.0-alpha); 
 }

 int Ndf=0;
 DF *df[10];

 int k=0;
 int Done=0;
 double kappa=dNaN;
 int OK, empty, kap_on;
 while (!Done)
 {
  df[Ndf]=new Std_DF(ne, ParmIn, k, &OK, &empty, &kap_on, &Done);
  if (!OK) 
  {
   IDLmsg("Invalid analytical distribution function parameters!");
   res=1;
  }

  if (OK && !empty) 
  {
   Ndf++;
   if (kap_on) kappa=ParmIn[4];
  }
  else delete df[Ndf];
  k++;
 }

 if (D1(ParmIn, InSize, 31))
 {
  df[Ndf]=new Arr_DF(ParmIn, InSize, NE, Nmu, E_arr, mu_arr, f3D, &OK, &empty);
  
  if (!OK)
  {
   IDLmsg("Invalid numerical distribution function parameters!");
   res=2;
  }

  if (OK && !empty) Ndf++;
  else delete df[Ndf];
 }
 
 df[Ndf]=0;
 
 if (!res)
 {
  double nb=0.0; //additional energetic electron density
  for (int i=0; i<Ndf; i++) nb+=df[i]->nb;

  *ne_total=ne+nb;
  double nu_p=e*sqrt(*ne_total/m/M_PI); 
  double nu_B=e*ParmIn[13]/m/c/(2.0*M_PI); 
  double theta=ParmIn[14]*M_PI/180;
  double nu_cr=nu_B*ParmIn[25];
  double nu_cr_WH=nu_B*ParmIn[26];
  int Npoints=(int)ParmIn[5];
  int m_on=(int)ParmIn[27];
  int Q_on=(int)ParmIn[28];

  for (int i=0; i<Nnu; i++) jX[i]=kX[i]=jO[i]=kO[i]=0;

  if (Ndf)
  {
   Find_jk_GS(df, nu, Nnu, -1, theta, nu_p, nu_B, nu_cr, nu_cr_WH, Npoints, Q_on, m_on, jX, kX);
   Find_jk_GS(df, nu, Nnu,  1, theta, nu_p, nu_B, nu_cr, nu_cr_WH, Npoints, Q_on, m_on, jO, kO);
  }

  if (T0>0)
  {
   Add_jk_FF(ne, nH, T0, kappa, nu, Nnu, -1, theta, nu_p, nu_B, jX, kX);
   Add_jk_FF(ne, nH, T0, kappa, nu, Nnu,  1, theta, nu_p, nu_B, jO, kO);
  }
 }

 for (int i=0; i<Ndf; i++) delete(df[i]);

 return res;
}

void RadiationTransfer(double nu, int Nsteps, double *dz, double *n0, double *B, double *theta, 
	                   double *jX, double *jO, double *kX, double *kO, 
					   double *Lw, double *Rw, double *Ls, double *Rs, double *Le, double *Re)
{
 for (int i=0; i<Nsteps; i++)
 {
  double tau=-kO[i]*dz[i];
  double eO=(tau<700) ? exp(tau) : 0.0; 
  double dIO=(kO[i]==0.0 || tau>700) ? 0.0 : jO[i]/kO[i]*((1.0-eO) ? 1.0-eO : -tau);
  tau=-kX[i]*dz[i];
  double eX=(tau<700) ? exp(tau) : 0.0;
  double dIX=(kX[i]==0.0 || tau>700) ? 0.0 : jX[i]/kX[i]*((1.0-eX) ? 1.0-eX : -tau);

  if (i>0) if ((theta[i]>(M_PI/2)) ^ (theta[i-1]>(M_PI/2)))
  {
   double a=*Lw;
   *Lw=*Rw;
   *Rw=a;
   
   double B_avg=(B[i]+B[i-1])/2;
   double n0_avg=(n0[i]+n0[i-1])/2;
   double da_dz=fabs(theta[i]-theta[i-1])/(dz[i]+dz[i-1])*2;

   double QT=e*e*e*e*e/(32*M_PI*M_PI*m*m*m*m*c*c*c*c)*n0_avg*sqr(B_avg)*B_avg/sqr(sqr(nu))/da_dz;
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

int MWTransfer(int Nsteps, int NE, int Nmu, int InSize, double *Parms, double *E_arr, double *mu_arr, double *f_arr, double *RL)
{
 #define OutSize 7

 int res=0;

 int Nnu=(int)D2(Parms, InSize, 18, 0);
 double *nu=(double*)malloc(sizeof(double)*Nnu);
 if (D2(Parms, InSize, 15, 0)>0)
 {
  nu[0]=D2(Parms, InSize, 15, 0);
  double dnu=D2(Parms, InSize, 16, 0);
  dnu=pow(10.0, dnu);
  for (int i=1; i<Nnu; i++) nu[i]=nu[i-1]*dnu;
 }
 else for (int i=0; i<Nnu; i++) nu[i]=D2(RL, OutSize, 0, i)*1e9;

 double **jX, **jO, **kX, **kO, *dz, *ne_total, *B, *theta;
 jX=(double**)malloc(sizeof(double*)*Nsteps);
 jO=(double**)malloc(sizeof(double*)*Nsteps);
 kX=(double**)malloc(sizeof(double*)*Nsteps);
 kO=(double**)malloc(sizeof(double*)*Nsteps);
 dz=(double*)malloc(sizeof(double)*Nsteps);
 ne_total=(double*)malloc(sizeof(double)*Nsteps);
 B=(double*)malloc(sizeof(double)*Nsteps);
 theta=(double*)malloc(sizeof(double)*Nsteps);
 for (int i=0; i<Nsteps; i++)
 {
  jX[i]=(double*)malloc(sizeof(double)*Nnu);
  jO[i]=(double*)malloc(sizeof(double)*Nnu);
  kX[i]=(double*)malloc(sizeof(double)*Nnu);
  kO[i]=(double*)malloc(sizeof(double)*Nnu);
 }

 int err=0;
 for (int i=0; i<Nsteps && !err; i++)
 {
  dz[i]=D2(Parms, InSize, 1, i);
  B[i]=D2(Parms, InSize, 13, i);
  theta[i]=D2(Parms, InSize, 14, i)*M_PI/180;

  IDL_DFarr *f3D=new IDL_DFarr(Nsteps, NE, Nmu, i, f_arr);
  err=FindLocalJK(nu, Nnu, Parms+i*InSize, InSize, 
	              NE, Nmu, E_arr, mu_arr, f3D, 
	              jX[i], jO[i], kX[i], kO[i], ne_total+i);
  delete f3D;
 }

 if (err)
 {
  IDLmsg("MWTransfer error: gyrosynchrotron calculation error.");
  res=err;
 }
 else
 {
  double Lw, Rw, Ls, Rs, Le, Re;

  double *jXl=(double*)malloc(sizeof(double)*Nsteps);
  double *jOl=(double*)malloc(sizeof(double)*Nsteps);
  double *kXl=(double*)malloc(sizeof(double)*Nsteps);
  double *kOl=(double*)malloc(sizeof(double)*Nsteps);

  double sang=D2(Parms, InSize, 0, 0)/(sqr(AU)*sfu);

  for (int i=0; i<Nnu; i++) 
  {
   Lw=D2(RL, OutSize, 1, i)/sang;
   Rw=D2(RL, OutSize, 2, i)/sang;
   Ls=D2(RL, OutSize, 3, i)/sang;
   Rs=D2(RL, OutSize, 4, i)/sang;
   Le=D2(RL, OutSize, 5, i)/sang;
   Re=D2(RL, OutSize, 6, i)/sang;

   for (int j=0; j<Nsteps; j++)
   {
	jXl[j]=jX[j][i];
	jOl[j]=jO[j][i];
	kXl[j]=kX[j][i];
	kOl[j]=kO[j][i];
   }

   RadiationTransfer(nu[i], Nsteps, dz, ne_total, B, theta, jXl, jOl, kXl, kOl, &Lw, &Rw, &Ls, &Rs, &Le, &Re);

   D2(RL, OutSize, 0, i)=nu[i]/1e9;
   D2(RL, OutSize, 1, i)=Lw*sang;
   D2(RL, OutSize, 2, i)=Rw*sang;
   D2(RL, OutSize, 3, i)=Ls*sang;
   D2(RL, OutSize, 4, i)=Rs*sang;
   D2(RL, OutSize, 5, i)=Le*sang;
   D2(RL, OutSize, 6, i)=Re*sang;
  }

  free(jXl);
  free(jOl);
  free(kXl);
  free(kOl);
 }

 for (int i=0; i<Nsteps; i++)
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