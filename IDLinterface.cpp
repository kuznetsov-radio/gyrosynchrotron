#include <stdlib.h>
#include <memory.h>
#include "MWmain.h"
#include "Plasma.h"
#include "IDLinterface.h"
#include "Messages.h"
#ifndef LINUX
#include <ppl.h>
#else
#include <omp.h>
#endif

#ifndef LINUX
extern "C" __declspec(dllexport) int GET_MW_main(int argc, void **argv)
#else
extern "C" double GET_MW_main(int argc, void **argv)
#endif
{
 int res=0;

 if (argc<7) res=-1;
 else
 {
  int *Lparms=(int*)argv[0];
  double *Rparms=(double*)argv[1];
  double *Parms=(double*)argv[2];
  double *E_arr=(double*)argv[3];
  double *mu_arr=(double*)argv[4];
  double *f_arr=(double*)argv[5];
  double *RL=(double*)argv[6];

  res=MW_Transfer(Lparms, Rparms, Parms, E_arr, mu_arr, f_arr, RL);
 }

 return res;
}

#ifndef LINUX
extern "C" __declspec(dllexport) int GET_MW(int argc, void **argv)
#else
extern "C" double GET_MW(int argc, void **argv)
#endif
{
 int res=0;

 if (argc<7) res=-1;
 else
 {
  int *Lparms=(int*)argv[0];
  double *Rparms=(double*)argv[1];
  double *Parms1=(double*)argv[2];
  double *E_arr=(double*)argv[3];
  double *mu_arr=(double*)argv[4];
  double *f_arr=(double*)argv[5];
  double *RL=(double*)argv[6];

  int Nz=Lparms[i_Nz];
  double *Parms=(double*)malloc(sizeof(double)*Nz*InSize);
  memcpy(Parms, Parms1, sizeof(double)*Nz*InSize);

  for (int i=0; i<Nz; i++)
  {
   double *p=Parms+i*InSize;

   if (p[i_T0]<1e5 && p[i_np]==0 && p[i_nH]==0)
   {
	double ne, nH, nHe;

	FindIonizationsSolar(p[i_n0], p[i_T0], &ne, &nH, &nHe);

	p[i_n0]=ne;
	p[i_nH]=nH;
	p[i_nHe]=nHe;
   }
  }

  res=MW_Transfer(Lparms, Rparms, Parms, E_arr, mu_arr, f_arr, RL);

  free(Parms);
 }

 return res;
}

#ifndef LINUX
extern "C" __declspec(dllexport) int GET_MW_SLICE(int argc, void **argv)
#else
extern "C" double GET_MW_SLICE(int argc, void **argv)
#endif
{
 int res=0;

 if (argc<7) res=-1;
 else
 {
  int *Lparms_M=(int*)argv[0];
  double *Rparms_M=(double*)argv[1];
  double *Parms_M=(double*)argv[2];
  double *E_arr=(double*)argv[3];
  double *mu_arr=(double*)argv[4];
  double *f_arr_M=(double*)argv[5];
  double *RL_M=(double*)argv[6];

  int Npix=Lparms_M[i_Npix];
  int Nz=Lparms_M[i_Nz+1];
  int Nnu=Lparms_M[i_Nnu+1];
  int NE=Lparms_M[i_NE+1];
  int Nmu=Lparms_M[i_Nmu+1];

  int *res_M=(int*)malloc(sizeof(int)*Npix);

  #ifndef LINUX

  concurrency::parallel_for(0, Npix, [&](int pix)
  {
   void *ARGV[7];
   ARGV[0]=(void*)(Lparms_M+1);
   ARGV[1]=(void*)(Rparms_M+pix*RpSize);
   ARGV[2]=(void*)(Parms_M+pix*Nz*InSize);
   ARGV[3]=(void*)E_arr;
   ARGV[4]=(void*)mu_arr;
   ARGV[5]=(void*)(f_arr_M+pix*Nz*NE*Nmu);
   ARGV[6]=(void*)(RL_M+pix*Nnu*OutSize);

   res_M[pix]=GET_MW(7, ARGV);
  });

  #else

  #pragma omp parallel for
  for(int pix=0; pix<Npix; pix++)
  {
   void *ARGV[7];
   ARGV[0]=(void*)(Lparms_M+1);
   ARGV[1]=(void*)(Rparms_M+pix*RpSize);
   ARGV[2]=(void*)(Parms_M+pix*Nz*InSize);
   ARGV[3]=(void*)E_arr;
   ARGV[4]=(void*)mu_arr;
   ARGV[5]=(void*)(f_arr_M+pix*Nz*NE*Nmu);
   ARGV[6]=(void*)(RL_M+pix*Nnu*OutSize);

   res_M[pix]=GET_MW(7, ARGV);
  }

  #endif

  for (int i=0; i<Npix; i++) res=(res!=0) || (res_M[i]!=0);

  free(res_M);
 }

 return res;
}

