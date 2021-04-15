#pragma once
#include "DF.h"

#define ERR_i 1e-5
#define ERR_s 1e-5
#define S_MIN 1
#define S_MAX 1000000
#define T_MAX 3600.0

void Find_jk_GS(DF **df, double *nu, int Nnu, int mode, double theta, double nu_p, double nu_B, 
	            double nu_cr, double nu_cr_WH, int Npoints, int Q_on, int m_on, double *j, double *k);