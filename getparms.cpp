#include <stdio.h>

const char* arr[]={
 "         dS;   0.180E+19      ;cm^2;                       Source/pixel Area",
 "         dR;   0.600E+09        ;cm;                      Source/voxel Depth",
 "        T_0;   0.200E+08         ;K;                      Plasma Temperature",
 "        eps;   0.500E-01      ;none;         Matching parm. for TNT distr-ns",
 "      kappa;    4.00          ;none;                  Index of Kappa-distr-n",
 "          N;          16      ;none;             Number of integration nodes",
 "       Emin;   0.100           ;MeV;                       Low energy cutoff",
 "       Emax;    10.0           ;MeV;                      High energy cutoff",
 "    E_break;    1.00           ;MeV;                         Break E for DPL",
 "     delta1;    4.00          ;none;                    (LE) Power-Law index",
 "     delta2;    6.00          ;none;            (HE) Power-Law index for DPL",
 "        n_0;   0.500E+10   ;cm^{-3};                       Thermal e density",
 "        n_b;   0.300E+08   ;cm^{-3};                    Nonthermal e density",
 "          B;    200.             ;G;                          Magnetic field",
 "      theta;    35.0       ;degrees;                           Viewing angle",
 "      f_min;   0.100E+10        ;Hz;        Starting freq. to calc. spectrum",
 "         df;   0.02        ;Log(Hz);           Logarithmic step in frequency",
 "     Dist_E;           3      ;none;        Type of distribution over energy",
 "     N_freq;         100      ;none;                   Number of frequencies",
 "   Dist_Ang;           1      ;none;            Type of angular distribution",
 "    theta_C;    60.0       ;degrees;                      Loss-cone boundary",
 "    theta_b;    90.0       ;degrees;                 Angle of beam direction",
 "        dMu;   0.100          ;none;         dMu for gau/exp/SuGau loss-cone",
 "        a_4;    10.0          ;none;       Coeff for a*(Mu-xMu0)^4 for SuGau",
 "   reserved;           0      ;none;                                reserved",
 "     f^C_cr;           0      ;f_ce;          Hybrid code boundary frequency",
 "    f^WH_cr;           0      ;f_ce;            Exact J_s boundary frequency",
 " Match. key;           1       ;1/0;                  Renormalization on/off",
 " Q-opt. key;           2     ;2/1/0;             Q-optimization on/appr./off",
 "        n_H;           0   ;cm^{-3};          Neutral hydrogen concentration",
 "        n_p;           0   ;cm^{-3};                   Protons concentration",
 "  Array key;           0       ;1/0;                            Array on/off",
 "     mu key;           2     ;2/1/0;      Exact (2) / Isotropic (1) / PK (0)",
 "    log key;           0       ;1/0;       Equidistant / logarithmic spacing",
 "   freq key;           0       ;1/0;         User-defined / log-spaced array",
 "   reserved;           0      ;none;                                reserved",
 " interp key;           0       ;1/0;               LQ / spline interpolation"
};

#define Nstrings 37

#ifndef LINUX
extern "C" __declspec(dllexport) float GET_PARMS(int argc, void **argv)
#else
extern "C" float GET_PARMS(int argc, void **argv)
#endif
{
 FILE *F=fopen("Parms.txt", "w");
 if (F)
 {
  for (int i=0; i<Nstrings; i++) fprintf(F, "%s\n", arr[i]);
  fclose(F);
  return 0;
 }
 else return -1;
}