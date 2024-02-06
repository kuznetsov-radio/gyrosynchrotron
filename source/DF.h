#pragma once

class DF
{
 public:
 int N_intervals{};
 double E_x[10]{};
 int logscale[9]{};
 double nb{}, EPS_mu0{};
 int PK_on;
 virtual void Fp(double p, double p_z, double p_n,
	             double *f, double *df_dp, double *df_dalpha)=0;
 virtual void FE(double E, double mu, 
	             double *f, double *df_dE, double *df_dmu, 
				 double *g1, double *g2)=0;
 DF()
 {
  PK_on=0;
 }
 virtual ~DF() {};
};

#define PosDef 1 // whether to require n_b > 0 or not