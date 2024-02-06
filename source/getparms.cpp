#include <stdio.h>

const char* arr1[]={
 "      N_vox;       1   ;int;     data;                    Number of voxels",
 "     N_freq;     100   ;int;     user;               Number of frequencies",
 "        N_E;       0   ;int;     data;                  Number of energies",
 "       N_mu;       0   ;int;     data;              Number of pitch_angles",
 "          N;       0   ;int;     user;         Number of integration nodes",
 " Match. key;       0   ;int;     user;              Renormalization on/off",
 " Q-opt. key;       0   ;int;     user;               Q-optimization on/off",
 "    arr_key;       0   ;int;     user;    Array DF on/off (for all voxels)",
 "    log_key;       0   ;int;     data;        Logarithmic E spacing on/off",
 "     PK_key;       0   ;int;     user;  Isotropic & Petrosian-Klein on/off",
 " spline_key;       0   ;int;     user;      Spline/LQ interpolation switch"
};

#define N1 11

const char* arr1s[]={
 "      N_pix;       1   ;int;     data;                    Number of pixels",
 "      N_vox;       1   ;int;     data;                    Number of voxels",
 "     N_freq;     100   ;int;     user;               Number of frequencies",
 "        N_E;       0   ;int;     data;                  Number of energies",
 "       N_mu;       0   ;int;     data;              Number of pitch_angles",
 "          N;       0   ;int;     user;         Number of integration nodes",
 " Match. key;       0   ;int;     user;              Renormalization on/off",
 " Q-opt. key;       0   ;int;     user;               Q-optimization on/off",
 "    arr_key;       0   ;int;     user;    Array DF on/off (for all voxels)",
 "    log_key;       0   ;int;     data;        Logarithmic E spacing on/off",
 "     PK_key;       0   ;int;     user;  Isotropic & Petrosian-Klein on/off",
 " spline_key;       0   ;int;     user;      Spline/LQ interpolation switch"
};

#define N1s 12

const char* arr2[]={
 "         dS;   1E+18   ;cm^2;    data;                   Source/pixel area",
 "      f_min;   1E+09   ;Hz;      user;    Starting freq. to calc. spectrum",
 "         df;    0.02   ;Log(Hz); user;       Logarithmic step in frequency",
 "     f^C_cr;       0   ;f_ce;    user;      Hybrid code boundary frequency",
 "    f^WH_cr;       0   ;f_ce;    user;        Exact J_s boundary frequency"
};

#define N2 5

const char *arr2d[]={
 "   distance;       1   ;AU;      user;            Distance from the source"
};

#define N2d 1

const char* arr3[]={
 "         dR;   1E+09   ;cm;      data;                  Source/voxel depth",
 "        T_0;   1E+06   ;K;       data;                  Plasma temperature",
 "        n_0;   1E+09   ;cm^{-3}; data;       Electron or total gas density",
 "          B;   200.0   ;G;       data;                      Magnetic field",
 "      theta;    60.0   ;degrees; data;                       Viewing angle",
 "  mech_flag;       0   ;int;     user;             Emission mechanism flag",
 "     Dist_E;       3   ;int;     user;    Type of distribution over energy",
 "        n_b;   1E+07   ;cm^{-3}; data;                Nonthermal e density",
 "  eps_kappa;       4   ;none;    user;     Parameter for TNT/kappa distrn.",
 "       Emin;     0.1   ;MeV;     user;                   Low-energy cutoff",
 "       Emax;    10.0   ;MeV;     user;                  High-energy cutoff",
 "    E_break;     1.0   ;MeV;     user;                Break energy for DPL",
 "     delta1;     4.0   ;none;    user;                (LE) Power-Law index",
 "     delta2;     6.0   ;none;    user;        (HE) Power-Law index for DPL",
 "   Dist_Ang;       1   ;none;    user;        Type of angular distribution",
 "    theta_0;    60.0   ;degrees; data; Loss-cone boundary / beam direction",
 "        dMu;     0.1   ;none;    user;     dMu for gau/exp/SuGau loss-cone",
 "        a_4;    10.0   ;none;    user;   Coeff for a*(Mu-xMu0)^4 for SuGau",
 "        n_p;       0   ;cm^{-3}; data;                Proton concentration",
 "       n_HI;       0   ;cm^{-3}; data;      Neutral hydrogen concentration",
 "      n_HeI;       0   ;cm^{-3}; data;        Neutral helium concentration",
 "arr_key_loc;       0   ;0/1;     data;                     Array DF on/off",
 "  abund_key;       0   ;-1/0/1;  data;           Dulk/Coronal/Photospheric",
 "   reserved;       0   ;int;     data;                            reserved"
};

#define N3 24

void WriteParms(const char **arr, const char *fname, int N, int add)
{
 FILE *F=fopen(fname, add ? "a" : "w");
 if (F)
 {
  for (int i=0; i<N; i++) fprintf(F, "%s\n", arr[i]);
  fclose(F);
 }
}

#ifndef LINUX
extern "C" __declspec(dllexport) float GET_PARMS(int argc, void **argv)
#else
extern "C" float GET_PARMS(int argc, void **argv)
#endif
{
 WriteParms(arr1, "Long_input.txt",  N1, 0);
 WriteParms(arr2, "Real_input.txt",  N2, 0);
 WriteParms(arr3, "Parms_input.txt", N3, 0);
 return 0;
}

#ifndef LINUX
extern "C" __declspec(dllexport) float GET_PARMS_SLICE(int argc, void **argv)
#else
extern "C" float GET_PARMS_SLICE(int argc, void **argv)
#endif
{
 WriteParms(arr1s, "Long_input.txt",  N1s, 0);
 WriteParms(arr2,  "Real_input.txt",  N2, 0);
 WriteParms(arr3,  "Parms_input.txt", N3, 0);
 return 0;
}

#ifndef LINUX
extern "C" __declspec(dllexport) float GET_PARMS1(int argc, void **argv)
#else
extern "C" float GET_PARMS1(int argc, void **argv)
#endif
{
 WriteParms(arr1, "Long_input.txt",  N1, 0);
 WriteParms(arr2, "Real_input.txt",  N2, 0);
 WriteParms(arr2d, "Real_input.txt",  N2d, 1);
 WriteParms(arr3, "Parms_input.txt", N3, 0);
 return 0;
}

#ifndef LINUX
extern "C" __declspec(dllexport) float GET_PARMS1_SLICE(int argc, void **argv)
#else
extern "C" float GET_PARMS1_SLICE(int argc, void **argv)
#endif
{
 WriteParms(arr1s, "Long_input.txt",  N1s, 0);
 WriteParms(arr2,  "Real_input.txt",  N2, 0);
 WriteParms(arr2d, "Real_input.txt",  N2d, 1);
 WriteParms(arr3,  "Parms_input.txt", N3, 0);
 return 0;
}