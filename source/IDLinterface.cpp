#include <stdlib.h>
#include <memory.h>
#include "MWmain.h"
#include "Plasma.h"
#include "IDLinterface.h"
#include "Messages.h"
#ifndef LINUX
#include <ppl.h>
#include <concrtrm.h>
#else
#include <omp.h>
#endif

#ifndef LINUX
extern "C" __declspec(dllexport) int GET_MW1_main(int argc, void **argv)
#else
extern "C" double GET_MW1_main(int argc, void **argv)
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
extern "C" __declspec(dllexport) int GET_MW1(int argc, void **argv)
#else
extern "C" double GET_MW1(int argc, void **argv)
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
extern "C" __declspec(dllexport) int GET_MW1_SLICE(int argc, void **argv)
#else
extern "C" double GET_MW1_SLICE(int argc, void **argv)
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
   
   res_M[pix]=GET_MW1(7, ARGV);
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
   
   res_M[pix]=GET_MW1(7, ARGV);
  }

  #endif

  for (int i=0; i<Npix; i++) res=(res!=0) || (res_M[i]!=0);

  free(res_M);
 }

 return res;
}

#ifndef LINUX
extern "C" __declspec(dllexport) int GET_MW_main(int argc,void** argv)
#else
extern "C" double GET_MW_main(int argc,void** argv)
#endif
{
 if (argc<7) return -1;
 else
 {
  double Rparms[RpSize];
  void* ARGV[7];

  memcpy(Rparms, argv[1], sizeof(double)*(RpSize-1));
  memcpy(ARGV, argv, sizeof(void*)*7);
  Rparms[i_dSun]=1;
  ARGV[1]=(void*)Rparms;

  return GET_MW1_main(7, ARGV);
 }
}

#ifndef LINUX
extern "C" __declspec(dllexport) int GET_MW(int argc, void** argv)
#else
extern "C" double GET_MW(int argc, void** argv)
#endif
{
 if (argc<7) return -1;
 else
 {
  double Rparms[RpSize];
  void* ARGV[7];

  memcpy(Rparms, argv[1], sizeof(double)*(RpSize-1));
  memcpy(ARGV, argv, sizeof(void*)*7);
  Rparms[i_dSun]=1;
  ARGV[1]=(void*)Rparms;

  return GET_MW1(7, ARGV);
 }
}

#ifndef LINUX
extern "C" __declspec(dllexport) int GET_MW_SLICE(int argc, void** argv)
#else
extern "C" double GET_MW_SLICE(int argc, void** argv)
#endif
{
 if (argc<7) return -1;
 else
 {
  int *Lparms_M=(int*)argv[0];
  int Npix=Lparms_M[i_Npix];

  double *Rparms_M=(double*)malloc(sizeof(double)*Npix*RpSize);
  for (int i=0; i<Npix; i++)
  {
   double *p=((double*)argv[1])+i*(RpSize-1);
   double *p1=Rparms_M+i*RpSize;
   memcpy(p1, p, sizeof(double)*(RpSize-1));
   p1[i_dSun]=1;
  }

  void *ARGV[7];
  memcpy(ARGV, argv, sizeof(void*)*7);
  ARGV[1]=(void*)Rparms_M;

  int res=GET_MW1_SLICE(7, ARGV);

  free(Rparms_M);

  return res;
 }
}

#ifndef LINUX
extern "C" __declspec(dllexport) int GET_MW_SLICE_M(int argc, void** argv)
#else
extern "C" double GET_MW_SLICE_M(int argc, void** argv)
#endif
{
 if (argc<7) return -1;
 else
 {
  int *Lparms_M=(int*)argv[0];
  int Nthreads=Lparms_M[12];

  #ifndef LINUX
  int NtMax=concurrency::GetProcessorCount();
  #else
  int NtMax=omp_get_max_threads();
  #endif

  if (Nthreads>NtMax) Nthreads=NtMax;

  if (Nthreads>0)
  {
   #ifndef LINUX
   if (concurrency::CurrentScheduler::Id()!=-1) concurrency::CurrentScheduler::Detach();
   concurrency::SchedulerPolicy policy;
   policy.SetConcurrencyLimits(Nthreads, Nthreads);
   concurrency::CurrentScheduler::Create(policy);
   #else
   omp_set_num_threads(Nthreads);
   #endif
  }

  return GET_MW_SLICE(7, argv);
 }
}