class EmWave
{
 public:
 int Valid, sigma;
 double nu, nu_p, nu_B, theta;
 double ct, st;
 double y, N, N_z, T, L, Zfactor;
 EmWave(double nu, double theta, int sigma, double nu_p, double nu_B);
};

#define cst_min 1e-3
#define S_MIN 1
#define S_MAX 1000000
#define T_MAX 3600.0
#define ERR_i 1e-5
#define ERR_s 1e-5

void GS_jk_mDF(EmWave *w, DF **df, int ExactBessel, double *j, double *k);
void GS_jk_approx_mDF(EmWave *w, DF **df, int Npoints, int Q_on, double *j, double *k);
void FF_jk(EmWave *w, double ne, double nH, double T0, double kappa, double *j, double *k);