#include "IDLinterface.h"
#include "Messages.h"

#ifndef LINUX
extern "C" __declspec(dllexport) int pyGET_MW(int *Lparms, double *Rparms, double *Parms,
                                              double *E_arr, double *mu_arr, double *f_arr, double *RL)
#else
extern "C" double pyGET_MW(int *Lparms, double *Rparms, double *Parms,
                           double *E_arr, double *mu_arr, double *f_arr, double *RL)
#endif
{
 void *ARGV[7];
 ARGV[0]=(void*)Lparms;
 ARGV[1]=(void*)Rparms;
 ARGV[2]=(void*)Parms;
 ARGV[3]=(void*)E_arr;
 ARGV[4]=(void*)mu_arr;
 ARGV[5]=(void*)f_arr;
 ARGV[6]=(void*)RL; 

 return GET_MW(7, ARGV);
}

#ifndef LINUX
extern "C" __declspec(dllexport) int pyGET_MW_SLICE(int *Lparms_M, double *Rparms_M, double *Parms_M,
                                                    double *E_arr, double *mu_arr, double *f_arr_M, 
                                                    double *RL_M)
#else
extern "C" double pyGET_MW_SLICE(int *Lparms_M, double *Rparms_M, double *Parms_M,
                                 double *E_arr, double *mu_arr, double *f_arr_M, double *RL_M)
#endif
{
 void *ARGV[7];
 ARGV[0]=(void*)Lparms_M;
 ARGV[1]=(void*)Rparms_M;
 ARGV[2]=(void*)Parms_M;
 ARGV[3]=(void*)E_arr;
 ARGV[4]=(void*)mu_arr;
 ARGV[5]=(void*)f_arr_M;
 ARGV[6]=(void*)RL_M;

 return GET_MW_SLICE(7, ARGV);
}