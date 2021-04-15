#pragma once
#include "DF.h"

#define FFF 0 //free-free only
#define FFF1 1 //same as free-free only
#define THM 2 //relativistic thermal distribution
#define PLW 3 //power-law distribution on energy
#define DPL 4 //double power-law distribution on energy
#define TNT 5 //thermal + nonthermal on energy
#define KAP 6 //kappa-distribution
#define PLP 7 //power-law distribution on impulse module
#define PLG 8 //power-law distribution on relativistic factor
#define TNP 9 //thermal + nonthermal on impulse module
#define TNG 10 //thermal + nonthermal on relativistic factor
#define TPL 11 //isotropic maxwellian + power-law on energy
#define TDP 12 //isotropic maxwellian + double power law

#define ISO 0 //isotropic distribution
#define ISO1 1 //same as isotropic
#define ELC 2 //exponential (on pitch-angle cosine) distribution
#define GAU 3 //gaussian (on pitch-angle cosine) distribution
#define GAB 4 //directed gaussian beam
#define SGA 5 //directed gaussian beam with 4th power term

class Std_DF_energy
{
 public:
 double nb;
 int N_intervals;
 double E_x[10];
 int logscale[9];
 virtual void FE(double E, double *f, double *f_E)=0;
 virtual void Fp(double p, double *f, double *f_p)=0;
 virtual ~Std_DF_energy() {};
};

class Std_DF_angle
{
 public:
 double EPS_mu0;
 virtual void Falpha(double mu, double sa, double *f, double *f_alpha)=0;
 virtual void Fmu(double mu, double *f, double *f_mu, double *g1, double *g2)=0;
 virtual double g1short(double mu)=0;
 Std_DF_angle()
 {
  EPS_mu0=1e-3;
 }
 virtual ~Std_DF_angle() {};
};

class Std_DF : public DF
{
 Std_DF_energy *F1;
 Std_DF_angle *F2;
 public:
 void Fp(double p, double p_z, double p_n, 
	     double *f, double *df_dp, double *df_dalpha);
 void FE(double E, double mu, 
	     double *f, double *df_dE, double *df_dmu, 
		 double *g1, double *g2);
 Std_DF(double *Parms, int k, int *OK, int *empty, int *kap_on, int *Done);
 ~Std_DF();
};