#pragma once
#include "DF.h"
#include "ExtMath.h"

class Arr_DF : public DF
{
 int NE{}, Nmu{};
 double *E_arr, *mu_arr, *f_avg;
 double **f_arr;
 int log_on{}, iso_on{}, LQ_on{};
 Spline *S1;
 Spline2D *S2;
 double IntegratedF();
 public:
 void Fp(double p, double p_z, double p_n, 
	     double *f, double *df_dp, double *df_dalpha);
 void FE(double E, double mu, 
	     double *f, double *df_dE, double *df_dmu, 
		 double *g1, double *g2);
 Arr_DF(int *Lparms, double *E, double *mu, double *f, int *OK, int *empty);
 ~Arr_DF();
};