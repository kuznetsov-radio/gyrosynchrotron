#pragma once

#define me 9.1093837015e-28
#define c 2.99792458e10
#define mc (me*c)
#define mc2 (me*c*c)
#define e (1.602176634e-19*c/10)
#define hPl 6.62607015e-27 //Planck constant h (without bar)
#define em_alpha (1.0/137.035999084) 
#define eV 1.602176634e-12
#define kB 1.380649e-16
#define AU 1.495978707e13 
#define sfu 1e-19
#define ieH 2.1798718e-11 //hydrogen ionization energy (theoretical), erg
#define ieHe12 3.9393356e-11 //helium 1st ionization energy, erg
#define ieHe2 8.71830663945283224e-11 //helium 2nd ionization energy, erg

#define cst_min 1e-3

class EmWave
{
 public:
 int Valid, sigma;
 double nu, nu_p, nu_B, theta;
 double ct{}, st{};
 double y{}, N{}, N_z{}, T{}, L{}, Zfactor{};
 EmWave(double nu, double theta, int sigma, double nu_p, double nu_B, int LT_on, int Zh_on);
};

void FindIonizationsSolar(double n0, double T0, double *n_e, double *n_H, double *n_He);