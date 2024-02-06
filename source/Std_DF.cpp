#include <math.h>
#include <float.h>
#include <memory.h>
#include "Plasma.h"
#include "ExtMath.h"
#include "Std_DF.h"
#include "IDLinterface.h"

class THMdf : public Std_DF_energy //2
{
 double A, theta;
 public:
 void FE(double E, double *f, double *f_E);
 void Fp(double p, double *f, double *f_p);
 THMdf(double *Parms, int *OK, int *empty);
 THMdf(double *Parms, int *OK, int *empty, double Emax);
};

void THMdf :: FE(double E, double *f, double *f_E)
{
 double G=E/mc2+1.0;
 double p=mc*sqrt(sqr(G)-1.0);

 double fp, dfp_dp;
 Fp(p, &fp, &dfp_dp);

 *f=fp*p*me*G;
 *f_E=sqr(me*G)*(dfp_dp+fp/p)+fp*p/(c*c);
}

void THMdf :: Fp(double p, double *f, double *f_p)
{
 double G=sqrt(1.0+sqr(p/mc));
 *f=A/(mc*mc*mc)*exp((1.0-G)/theta);
 *f_p=-(*f)*p/G/theta/(mc*mc);
}

THMdf :: THMdf(double *Parms, int *OK, int *empty)
{
 double T0=fabs(Parms[i_T0]);

 theta=kB*T0/mc2;
 A=Parms[i_n0]/(2.0*M_PI*theta*ExpBesselK(2, 1.0/theta));

 N_intervals=3;
 E_x[0]=0.0;
 E_x[1]=kB*T0;
 E_x[2]=kB*T0*3;
 E_x[3]=kB*T0*710;
 logscale[0]=0;
 logscale[1]=0;
 logscale[2]=1;

 *OK=finite(A)!=0 && A>=0;
 *empty=(A==0.0);
 nb=0.0;
}

THMdf :: THMdf(double *Parms, int *OK, int *empty, double Emax)
{
 double T0=fabs(Parms[i_T0]);

 theta=kB*T0/mc2;
 A=Parms[i_n0]/(2.0*M_PI*theta*ExpBesselK(2, 1.0/theta));

 N_intervals=3;
 E_x[0]=0.0;
 E_x[1]=kB*T0;
 E_x[2]=kB*T0*3;
 E_x[3]=kB*T0*710;
 logscale[0]=0;
 logscale[1]=0;
 logscale[2]=1;

 if (Emax<=E_x[1])
 {
  N_intervals=1;
  E_x[1]=Emax;
 }
 else if (Emax<=E_x[2])
 {
  N_intervals=2;
  E_x[2]=Emax;
 }
 else if (Emax<=E_x[3]) E_x[3]=Emax;

 *OK=finite(A)!=0 && A>=0;
 *empty=(A==0.0);
 nb=0.0;
}

class PLWdf : public Std_DF_energy //3
{
 double A, delta;
 public:
 void FE(double E, double *f, double *f_E);
 void Fp(double p, double *f, double *f_p);
 PLWdf(double *Parms, int *OK, int *empty);
};

void PLWdf :: FE(double E, double *f, double *f_E)
{
 *f=A*pow(E, -delta);
 *f_E=(-delta/E)*(*f);
}

void PLWdf :: Fp(double p, double *f, double *f_p)
{
 double G=sqrt(1.0+sqr(p/mc));
 double E=mc2*(G-1.0);

 double fE, dfE_dE;
 FE(E, &fE, &dfE_dE);

 *f=fE/(p*me*G);
 *f_p=(dfE_dE-fE*G*me/sqr(p)*(1.0+sqr(p/G/mc)))/(sqr(me*G));
}

PLWdf :: PLWdf(double *Parms, int *OK, int *empty)
{
 double E1=Parms[i_Emin]*eV*1e6;
 double E2=Parms[i_Emax]*eV*1e6;
 delta=Parms[i_delta1];
 nb=Parms[i_nb];

 A=nb/(2.0*M_PI)*(delta-1.0)/(pow(E1, 1.0-delta)-pow(E2, 1.0-delta));

 N_intervals=1;
 E_x[0]=E1;
 E_x[1]=E2;
 logscale[0]=1;

 *OK=finite(A)!=0 && (PosDef ? A>=0.0 : 1) && E2>E1;
 *empty=(nb==0.0);
}

class DPLdf : public Std_DF_energy //4
{
 double Ebr, A1, A2, delta1, delta2;
 public:
 void FE(double E, double *f, double *f_E);
 void Fp(double p, double *f, double *f_p);
 DPLdf(double *Parms, int *OK, int *empty);
};

void DPLdf :: FE(double E, double *f, double *f_E)
{
 if (E<Ebr)
 {
  *f=A1*pow(E, -delta1);
  *f_E=(-delta1/E)*(*f);
 }
 else
 {
  *f=A2*pow(E, -delta2);
  *f_E=(-delta2/E)*(*f);
 }
}

void DPLdf :: Fp(double p, double *f, double *f_p)
{
 double G=sqrt(1.0+sqr(p/mc));
 double E=mc2*(G-1.0);

 double fE, dfE_dE;
 FE(E, &fE, &dfE_dE);

 *f=fE/(p*me*G);
 *f_p=(dfE_dE-fE*G*me/sqr(p)*(1.0+sqr(p/G/mc)))/(sqr(me*G));
}

DPLdf :: DPLdf(double *Parms, int *OK, int *empty)
{
 double E1=Parms[i_Emin]*eV*1e6;
 double E2=Parms[i_Emax]*eV*1e6;
 Ebr=Parms[i_Ebreak]*eV*1e6;
 delta1=Parms[i_delta1];
 delta2=Parms[i_delta2];
 nb=Parms[i_nb];

 A1=nb/(2.0*M_PI)/
	((pow(E1, 1.0-delta1)-pow(Ebr, 1.0-delta1))/(delta1-1.0)+
	 pow(Ebr, delta2-delta1)*(pow(Ebr, 1.0-delta2)-pow(E2, 1.0-delta2))/(delta2-1.0));
 A2=A1*pow(Ebr, delta2-delta1);

 N_intervals=2;
 E_x[0]=E1;
 E_x[1]=Ebr;
 E_x[2]=E2;
 logscale[0]=1;
 logscale[1]=1;

 *OK=finite(A1)!=0 && (PosDef ? A1>=0.0 : 1) && 
	 finite(A2)!=0 && (PosDef ? A2>=0.0 : 1) && 
	 E2>Ebr && Ebr>E1;
 *empty=(nb==0.0);
}

class KAPdf : public Std_DF_energy //6
{
 double A, kappa_p1, kappa_m32_theta;
 public:
 void FE(double E, double *f, double *f_E);
 void Fp(double p, double *f, double *f_p);
 KAPdf(double *Parms, int *OK, int *empty);
};

void KAPdf :: FE(double E, double *f, double *f_E)
{
 double G=E/mc2+1.0;
 double p=mc*sqrt(sqr(G)-1.0);

 double fp, dfp_dp;
 Fp(p, &fp, &dfp_dp);

 *f=fp*p*me*G;
 *f_E=sqr(me*G)*(dfp_dp+fp/p)+fp*p/(c*c);
}

void KAPdf :: Fp(double p, double *f, double *f_p)
{
 double G=sqrt(1.0+sqr(p/mc));
 double D=1.0+(G-1.0)/kappa_m32_theta;
 *f=A/(mc*mc*mc)*pow(D, -kappa_p1);
 *f_p=-kappa_p1*(*f)/D/kappa_m32_theta*p/G/(mc*mc);
}

class KappaIntegrand : public IntegrableFunction
{
 public:
 double kappa_m32_theta{}, kappa_p1{};
 double F(double G);
};

double KappaIntegrand :: F(double G)
{
 return G*sqrt(sqr(G)-1.0)*pow(1.0+(G-1.0)/kappa_m32_theta, -kappa_p1);
}

KAPdf :: KAPdf(double *Parms, int *OK, int *empty)
{
 double T0=fabs(Parms[i_T0]);
 double kappa=Parms[i_epskappa];
 double E_max=Parms[i_Emax]*eV*1e6;

 double theta=kB*T0/mc2;
 kappa_p1=kappa+1.0;
 kappa_m32_theta=(kappa-1.5)*theta;

 double G_max=E_max/mc2+1.0;
 KappaIntegrand ki;
 ki.kappa_m32_theta=kappa_m32_theta;
 ki.kappa_p1=kappa_p1;
 int err;
 A=Parms[i_n0]/(2.0*M_PI)/qromb(&ki, 1.0, G_max, 1e-6, &err);

 E_x[0]=0.0;
 if (E_max<=(kB*T0))
 {
  N_intervals=1;
  E_x[1]=E_max;
  logscale[0]=0;
 }
 else if (E_max<=(kB*T0*3))
 {
  N_intervals=2;
  E_x[1]=kB*T0;
  E_x[2]=E_max;
  logscale[0]=logscale[1]=0;
 }
 else
 {
  N_intervals=3;
  E_x[1]=kB*T0;
  E_x[2]=kB*T0*3;
  E_x[3]=E_max;
  logscale[0]=logscale[1]=0;
  logscale[2]=(E_x[3]/E_x[2]>5.0);
 }

 *OK=finite(A)!=0 && A>=0.0;
 *empty=(A==0.0);
 nb=0.0;
}

class PLPdf : public Std_DF_energy //7
{
 double A, delta;
 public:
 void FE(double E, double *f, double *f_E);
 void Fp(double p, double *f, double *f_p);
 PLPdf(double *Parms, int *OK, int *empty);
};

void PLPdf :: FE(double E, double *f, double *f_E)
{
 double G=E/mc2+1.0;
 double p=mc*sqrt(sqr(G)-1.0);

 double fp, dfp_dp;
 Fp(p, &fp, &dfp_dp);

 *f=fp*p*me*G;
 *f_E=sqr(me*G)*(dfp_dp+fp/p)+fp*p/(c*c);
}

void PLPdf :: Fp(double p, double *f, double *f_p)
{
 *f=A*pow(p, -delta);
 *f_p=-delta*(*f)/p;
}

PLPdf :: PLPdf(double *Parms, int *OK, int *empty)
{
 double E1=Parms[i_Emin]*eV*1e6;
 double E2=Parms[i_Emax]*eV*1e6;
 delta=Parms[i_delta1];
 nb=Parms[i_nb];

 double p1=mc*sqrt(sqr(E1/mc2+1.0)-1.0);
 double p2=mc*sqrt(sqr(E2/mc2+1.0)-1.0);

 A=nb/(2.0*M_PI)*(delta-3.0)/(pow(p1, 3.0-delta)-pow(p2, 3.0-delta));

 N_intervals=1;
 E_x[0]=E1;
 E_x[1]=E2;
 logscale[0]=1;

 *OK=finite(A)!=0 && (PosDef ? A>=0.0 : 1) && E2>E1;
 *empty=(nb==0.0);
}

class PLGdf : public Std_DF_energy //8
{
 double A, delta;
 public:
 void FG(double G, double *f, double *f_G);
 void FE(double E, double *f, double *f_E);
 void Fp(double p, double *f, double *f_p);
 PLGdf(double *Parms, int *OK, int *empty);
};

void PLGdf :: FG(double G, double *f, double *f_G)
{
 *f=A*pow(G, -delta);
 *f_G=-delta*(*f)/G;
}

void PLGdf :: FE(double E, double *f, double *f_E)
{
 double G=E/mc2+1.0;

 double fG, dfG_dG;
 FG(G, &fG, &dfG_dG);

 *f=fG/mc2;
 *f_E=dfG_dG/(mc2*mc2);
}

void PLGdf :: Fp(double p, double *f, double *f_p)
{
 double G=sqrt(1.0+sqr(p/mc));

 double fG, dfG_dG;
 FG(G, &fG, &dfG_dG);

 *f=fG/(p*G*mc*mc);
 *f_p=(dfG_dG-fG/G*(sqr(G*mc/p)+1.0))/sqr(G*mc*mc);
}

PLGdf :: PLGdf(double *Parms, int *OK, int *empty)
{
 double E1=Parms[i_Emin]*eV*1e6;
 double E2=Parms[i_Emax]*eV*1e6;
 delta=Parms[i_delta1];
 nb=Parms[i_nb];

 double G1=E1/mc2+1.0;
 double G2=E2/mc2+1.0;

 A=nb/(2.0*M_PI)*(delta-1.0)/(pow(G1, 1.0-delta)-pow(G2, 1.0-delta));

 N_intervals=1;
 E_x[0]=E1;
 E_x[1]=E2;
 logscale[0]=1;

 *OK=finite(A)!=0 && (PosDef ? A>=0.0 : 1) && E2>E1;
 *empty=(nb==0.0);
}

class ISOdf : public Std_DF_angle //0, 1
{
 public:
 void Falpha(double mu, double sa, double *f, double *f_alpha);
 void Fmu(double mu, double *f, double *f_mu, double *g1, double *g2);
 double g1short(double mu);
 ISOdf(int *OK);
};

void ISOdf :: Falpha(double mu, double sa, double *f, double *f_alpha)
{
 *f=0.5;
 *f_alpha=0.0;
}

void ISOdf :: Fmu(double mu, double *f, double *f_mu, double *g1, double *g2)
{
 *f=0.5;
 *f_mu=0.0;
 if (g1) *g1=0.0;
 if (g2) *g2=0.0;
}

double ISOdf :: g1short(double mu)
{
 return 0.0;
}

ISOdf :: ISOdf(int *OK)
{
 *OK=1;
}

class ELCdf : public Std_DF_angle //2
{
 double B, alpha_c, mu_c, dmu;
 public:
 void Falpha(double mu, double sa, double *f, double *f_alpha);
 void Fmu(double mu, double *f, double *f_mu, double *g1, double *g2);
 double g1short(double mu);
 ELCdf(double *Parms, int *OK);
};

void ELCdf :: Falpha(double mu, double sa, double *f, double *f_alpha)
{
 Fmu(mu, f, f_alpha, 0, 0);
 *f_alpha*=(-sa);
}

void ELCdf :: Fmu(double mu, double *f, double *f_mu, double *g1, double *g2)
{
 double amu=fabs(mu);
 double g1loc;
 if (amu<mu_c)
 {
  *f=B;
  *f_mu=0.0;
  if (g1) *g1=0.0;
  if (g2) *g2=0.0;
 }
 else
 {
  *f=B*exp(-(amu-mu_c)/dmu);
  g1loc=-sign(mu)/dmu;
  *f_mu=*f*g1loc;
  if (g1) *g1=g1loc;
  if (g2) *g2=1.0/sqr(dmu);
 }
}

double ELCdf :: g1short(double mu)
{
 double amu=fabs(mu);
 if (amu<mu_c) return 0.0;
 else return -sign(mu)/dmu;
}

ELCdf :: ELCdf(double *Parms, int *OK)
{
 alpha_c=Parms[i_alphaC]*M_PI/180;
 dmu=Parms[i_dmu];

 mu_c=fabs(cos(alpha_c));
 B=0.5/(mu_c+dmu-dmu*exp((mu_c-1.0)/dmu));

 EPS_mu0=min(EPS_mu0, dmu/30);
 *OK=finite(B)!=0 && B>0.0;
}

class GAUdf : public Std_DF_angle //3
{
 double B, alpha_c, mu_c, dmu;
 public:
 void Falpha(double mu, double sa, double *f, double *f_alpha);
 void Fmu(double mu, double *f, double *f_mu, double *g1, double *g2);
 double g1short(double mu);
 GAUdf(double *Parms, int *OK);
};

void GAUdf :: Falpha(double mu, double sa, double *f, double *f_alpha)
{
 Fmu(mu, f, f_alpha, 0, 0);
 *f_alpha*=(-sa);
}

void GAUdf :: Fmu(double mu, double *f, double *f_mu, double *g1, double *g2)
{
 double amu=fabs(mu);
 double g1loc;
 if (amu<mu_c)
 {
  *f=B;
  *f_mu=0.0;
  if (g1) *g1=0.0;
  if (g2) *g2=0.0;
 }
 else
 {
  *f=B*exp(-sqr((amu-mu_c)/dmu));
  g1loc=-2.0*(amu-mu_c)/sqr(dmu)*sign(mu);
  *f_mu=*f*g1loc;
  if (g1) *g1=g1loc;
  if (g2) *g2=4.0*sqr((amu-mu_c)/sqr(dmu))-2.0/sqr(dmu);
 }
}

double GAUdf :: g1short(double mu)
{
 double amu=fabs(mu);
 if (amu<mu_c) return 0.0;
 else return -2.0*(amu-mu_c)/sqr(dmu)*sign(mu);
}

GAUdf :: GAUdf(double *Parms, int *OK)
{
 alpha_c=Parms[i_alphaC]*M_PI/180;
 dmu=Parms[i_dmu];

 mu_c=fabs(cos(alpha_c));
 B=0.5/(mu_c+dmu*sqrt(M_PI)/2*Erf((1.0-mu_c)/dmu));

 EPS_mu0=min(EPS_mu0, sqr(dmu)/30);
 *OK=finite(B)!=0 && B>0.0;
}

class GABdf : public Std_DF_angle //4
{
 double B, alpha_c, mu_c, dmu;
 public:
 void Falpha(double mu, double sa, double *f, double *f_alpha);
 void Fmu(double mu, double *f, double *f_mu, double *g1, double *g2);
 double g1short(double mu);
 GABdf(double *Parms, int *OK);
};

void GABdf :: Falpha(double mu, double sa, double *f, double *f_alpha)
{
 Fmu(mu, f, f_alpha, 0, 0);
 *f_alpha*=(-sa);
}

void GABdf :: Fmu(double mu, double *f, double *f_mu, double *g1, double *g2)
{
 *f=B*exp(-sqr((mu-mu_c)/dmu));
 double g1loc=-2.0*(mu-mu_c)/sqr(dmu);
 *f_mu=*f*g1loc;
 if (g1) *g1=g1loc;
 if (g2) *g2=4.0*sqr((mu-mu_c)/sqr(dmu))-2.0/sqr(dmu);
}

double GABdf :: g1short(double mu)
{
 return -2.0*(mu-mu_c)/sqr(dmu);
}

GABdf :: GABdf(double *Parms, int *OK)
{
 alpha_c=Parms[i_alphaC]*M_PI/180;
 dmu=Parms[i_dmu];

 mu_c=cos(alpha_c);
 B=2.0/(sqrt(M_PI)*dmu)/(Erf((1.0-mu_c)/dmu)+Erf((1.0+mu_c)/dmu));

 EPS_mu0=min(EPS_mu0, sqr(dmu)/30);
 *OK=finite(B)!=0 && B>0.0;
}

class SGAdf : public Std_DF_angle //5
{
 double B, alpha_c, mu_c, dmu, a4;
 public:
 void Falpha(double mu, double sa, double *f, double *f_alpha);
 void Fmu(double mu, double *f, double *f_mu, double *g1, double *g2);
 double g1short(double mu);
 SGAdf(double *Parms, int *OK);
};

void SGAdf :: Falpha(double mu, double sa, double *f, double *f_alpha)
{
 Fmu(mu, f, f_alpha, 0, 0);
 *f_alpha*=(-sa);
}

void SGAdf :: Fmu(double mu, double *f, double *f_mu, double *g1, double *g2)
{
 double d2=sqr(mu-mu_c);
 *f=B*exp(-(d2+a4*sqr(d2))/sqr(dmu));
 double g1loc=-2.0*(mu-mu_c)*(1.0+2.0*a4*d2)/sqr(dmu);
 *f_mu=*f*g1loc;
 if (g1) *g1=g1loc;
 if (g2) *g2=2.0/sqr(sqr(dmu))*
	         (2.0*d2*(1.0-3.0*a4*sqr(dmu)+4.0*a4*d2+4.0*sqr(a4*d2))-sqr(dmu)); 
}

double SGAdf :: g1short(double mu)
{
 double d2=sqr(mu-mu_c);
 return -2.0*(mu-mu_c)*(1.0+2.0*a4*d2)/sqr(dmu);
}

class SGAIntegrand : public IntegrableFunction
{
 public:
 double mu_c{}, dmu{}, a4{};
 double F(double mu);
};

double SGAIntegrand :: F(double mu)
{
 double d2=sqr(mu-mu_c);
 return exp(-(d2+a4*sqr(d2))/sqr(dmu));
}

SGAdf :: SGAdf(double *Parms, int *OK)
{
 alpha_c=Parms[i_alphaC]*M_PI/180;
 dmu=Parms[i_dmu];
 a4=Parms[i_a4];

 mu_c=cos(alpha_c);
 SGAIntegrand sgi;
 sgi.mu_c=mu_c;
 sgi.dmu=dmu;
 sgi.a4=a4;
 int err;
 B=1.0/qromb(&sgi, -1, 1, 1e-10, &err);

 EPS_mu0=min(EPS_mu0, sqr(dmu)/30);
 *OK=finite(B)!=0 && B>0.0;
}

void Std_DF :: Fp(double p, double p_z, double p_n, double *f, double *df_dp, double *df_dalpha)
{
 double f1, f2, df1_dp, df2_dalpha;

 F1->Fp(p, &f1, &df1_dp);

 double mu=(p>0.0) ? p_z/p : 0.0;
 double sa=(p>0.0) ? p_n/p : 0.0;
 if (mu>1.0) mu=1.0;
 if (mu<(-1.0)) mu=-1.0;
 if (sa>1.0) sa=1.0;
 if (sa<(-1.0)) sa=-1.0;
  
 F2->Falpha(mu, sa, &f2, &df2_dalpha);

 *f=f1*f2; 
 *df_dp=df1_dp*f2;
 *df_dalpha=f1*df2_dalpha;
}

void Std_DF :: FE(double E, double mu, 
                  double *f, double *df_dE, double *df_dmu, 
		          double *g1, double *g2)
{
 if (!f) *g1=F2->g1short(mu); //only g1 is needed
 else //calculating all
 {
  double f1, f2, df1_dE, df2_dmu;

  F1->FE(E, &f1, &df1_dE);
  F2->Fmu(mu, &f2, &df2_dmu, g1, g2);

  *f=f1*f2;
  *df_dE=df1_dE*f2;
  *df_dmu=f1*df2_dmu;
 }
}

Std_DF :: Std_DF(double *Parms, int k, int *OK, int *empty, int *kap_on, int *Done)
{
 F1=0;
 F2=0;

 int E_id=(int)Parms[i_EId];
 int mu_id=(int)Parms[i_muId];

 *kap_on=0;

 switch (E_id)
 {
  case THM: F1=new THMdf(Parms, OK, empty);
	        *Done=1;
			break;

  case PLW: F1=new PLWdf(Parms, OK, empty);
	        *Done=1;
			break;

  case DPL: F1=new DPLdf(Parms, OK, empty);
	        *Done=1;
			break;

  case TNT: 
  case TNP:
  case TNG: {
	         double pcr=mc*sqrt((sqr(kB*Parms[i_T0]/mc2+1.0)-1.0)/Parms[i_epskappa]);
	         double Gcr=sqrt(1.0+sqr(pcr/mc));
	         double Ecr=mc2*(Gcr-1.0);
	  
	         if (!k) F1=new THMdf(Parms, OK, empty, Ecr);
			 else
			 {
	          double E_max=Parms[i_Emax]*eV*1e6;
			  double G_max=E_max/mc2+1.0;
              double p_max=mc*sqrt(sqr(G_max)-1.0);
			  double delta=Parms[i_delta1];

			  if (E_max>Ecr)
			  {
			   Std_DF_energy *thm=new THMdf(Parms, OK, empty);

			   if (*OK && !(*empty))
			   {
			    double fcr, tmp, Acr, nb;

				switch (E_id)
				{
				 case TNT: thm->FE(Ecr, &fcr, &tmp);
			               Acr=fcr*pow(Ecr, delta);
			               nb=Acr*2.0*M_PI/(delta-1.0)*(pow(Ecr, 1.0-delta)-pow(E_max, 1.0-delta));
						   break;
				 case TNP: thm->Fp(pcr, &fcr, &tmp);
					       Acr=fcr*pow(pcr, delta);
						   nb=Acr*2.0*M_PI/(delta-3.0)*(pow(pcr, 3.0-delta)-pow(p_max, 3.0-delta));
						   break;
				 case TNG: thm->FE(Ecr, &fcr, &tmp);
					       Acr=fcr*mc2*pow(Gcr, delta);
						   nb=Acr*2.0*M_PI/(delta-1.0)*(pow(Gcr, 1.0-delta)-pow(G_max, 1.0-delta));
				}

			    double ParmsLoc[InSize];
				memcpy(ParmsLoc, Parms, sizeof(double)*InSize);
			    ParmsLoc[i_Emin]=Ecr/eV/1e6;
                ParmsLoc[i_nb]=nb;

				switch (E_id)
				{
				 case TNT: F1=new PLWdf(ParmsLoc, OK, empty);
					       break;
                 case TNP: F1=new PLPdf(ParmsLoc, OK, empty);
					       break;
                 case TNG: F1=new PLGdf(ParmsLoc, OK, empty);
				}
			   }

			   delete thm;
			  }
			  else
			  {
			   *OK=1;
			   *empty=1;
			  }

			  *Done=1;
			 }
			}
			break;

  case KAP: F1=new KAPdf(Parms, OK, empty);
	        *Done=1;
			*kap_on=1;
			break;

  case PLP: F1=new PLPdf(Parms, OK, empty);
	        *Done=1;
			break;

  case PLG: F1=new PLGdf(Parms, OK, empty);
	        *Done=1;
			break;

  case TPL:
  case TDP: if (!k)
			{
			 F1=new THMdf(Parms, OK, empty);
			 mu_id=ISO;
			}
			else
			{
			 switch (E_id)
			 {
			  case TPL: F1=new PLWdf(Parms, OK, empty);
				        break;
			  case TDP: F1=new DPLdf(Parms, OK, empty);
			 }

			 *Done=1;
			}
			break;

  default: *OK=1;
	       *empty=1;
		   *Done=1;
 }

 if (*OK && !(*empty))
 {
  nb=F1->nb;
  N_intervals=F1->N_intervals;
  for (int i=0; i<N_intervals; i++) logscale[i]=F1->logscale[i];
  for (int i=0; i<=N_intervals; i++) E_x[i]=F1->E_x[i];

  switch (mu_id)
  {
   case ISO:
   case ISO1: F2=new ISOdf(OK);
	          break;

   case ELC: F2=new ELCdf(Parms, OK);
	         break;

   case GAU: F2=new GAUdf(Parms, OK);
	         break;

   case GAB: F2=new GABdf(Parms, OK);
	         break;

   case SGA: F2=new SGAdf(Parms, OK);
	         break;

   default: F2=new ISOdf(OK);
  }

  if (*OK) EPS_mu0=F2->EPS_mu0;
 }
}

Std_DF :: ~Std_DF()
{
 if (F1) delete F1;
 if (F2) delete F2;
}