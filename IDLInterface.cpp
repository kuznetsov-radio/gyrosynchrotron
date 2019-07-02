#include "Messages.h"
#include "MWmain.h"
#include "ExtMath.h"
#include "IDLInterface.h"
#ifndef LINUX
#include <ppl.h>
#else
#include <omp.h>
#endif

#ifndef LINUX
extern "C" __declspec(dllexport) double GET_MW(int argc, void **argv)
#else
extern "C" double GET_MW(int argc, void **argv)
#endif
{
 if (argc<6)
 {
  IDLmsg("GET_MW error: not enough parameters in the function call.");
  return -1;
 }

 int *Ndat=(int*)argv[0];
 double *Parms=(double*)argv[1];
 double *E_arr=(double*)argv[2];
 double *mu_arr=(double*)argv[3];
 double *f_arr=(double*)argv[4];
 double *RL=(double*)argv[5];

 if (Ndat[0]<=0)
 {
  IDLmsg("GET_MW error: number of voxels must be positive.");
  return -2;
 }

 if (Ndat[3]<29) //classical interface - ParmIn=fltarr(29)
 {
  IDLmsg("GET_MW error: not enough elements in the voxel parameters array.");
  return -2;
 }

 int res=MWTransfer(Ndat[0], Ndat[1], Ndat[2], Ndat[3], Parms, E_arr, mu_arr, f_arr, RL);

 return (double)res;
}

#ifndef LINUX
extern "C" __declspec(dllexport) double GET_MW_SLICE(int argc, void **argv)
#else
extern "C" double GET_MW_SLICE(int argc, void **argv)
#endif
{
 #define OutSize 7

 if (argc<6)
 {
  IDLmsg("GET_MW_SLICE error: not enough parameters in the function call.");
  return -1;
 }

 int *Ndat=(int*)argv[0];
 double *Parms=(double*)argv[1];
 double *E_arr=(double*)argv[2];
 double *mu_arr=(double*)argv[3];
 double *f_arr=(double*)argv[4];
 double *RL=(double*)argv[5];

 int Npix=Ndat[0];
 int Nsteps=Ndat[1];
 int NE=Ndat[2];
 int Nmu=Ndat[3];
 int InSize=Ndat[4];
 int Nnu=(int)D3(Parms, InSize, Nsteps, 18, 0, 0);

 if (Npix<=0)
 {
  IDLmsg("GET_MW_SLICE error: number of pixels must be positive.");
  return -2;
 }

 #ifndef LINUX
 concurrency::parallel_for(0, Npix, [&](int pix)
 {
  if (pix>0)
  {
	D3(Parms, InSize, Nsteps, 15, 0, pix)=D3(Parms, InSize, Nsteps, 15, 0, 0);
	D3(Parms, InSize, Nsteps, 16, 0, pix)=D3(Parms, InSize, Nsteps, 16, 0, 0);
	D3(Parms, InSize, Nsteps, 18, 0, pix)=D3(Parms, InSize, Nsteps, 18, 0, 0);
	for (int i=0; i<Nnu; i++) D3(RL, OutSize, Nnu, 0, i, pix)=D3(RL, OutSize, Nnu, 0, i, 0);
   }

   void *ARGV[6];
   ARGV[0]=(void*)(Ndat+1);
   ARGV[1]=(void*)(Parms+pix*Nsteps*InSize);
   ARGV[2]=(void*)(E_arr);
   ARGV[3]=(void*)(mu_arr);
   ARGV[4]=(void*)(f_arr+pix*Nsteps*NE*Nmu);
   ARGV[5]=(void*)(RL+pix*Nnu*OutSize);

   GET_MW(6, ARGV);
 });
 #else
 #pragma omp parallel for
 for (int pix=0; pix<Npix; pix++)
 {
  if (pix>0)
  {
	D3(Parms, InSize, Nsteps, 15, 0, pix)=D3(Parms, InSize, Nsteps, 15, 0, 0);
	D3(Parms, InSize, Nsteps, 16, 0, pix)=D3(Parms, InSize, Nsteps, 16, 0, 0);
	D3(Parms, InSize, Nsteps, 18, 0, pix)=D3(Parms, InSize, Nsteps, 18, 0, 0);
	for (int i=0; i<Nnu; i++) D3(RL, OutSize, Nnu, 0, i, pix)=D3(RL, OutSize, Nnu, 0, i, 0);
   }

   void *ARGV[6];
   ARGV[0]=(void*)(Ndat+1);
   ARGV[1]=(void*)(Parms+pix*Nsteps*InSize);
   ARGV[2]=(void*)(E_arr);
   ARGV[3]=(void*)(mu_arr);
   ARGV[4]=(void*)(f_arr+pix*Nsteps*NE*Nmu);
   ARGV[5]=(void*)(RL+pix*Nnu*OutSize);

   GET_MW(6, ARGV);
 }
 #endif

 return 0;
}

IDL_DFarr :: IDL_DFarr(int _Nsteps, int _NE, int _Nmu, int _step, double *_f)
{
 Nsteps=_Nsteps;
 NE=_NE;
 Nmu=_Nmu;
 step=_step;
 f=_f;
}

double IDL_DFarr :: A(int i, int j)
{
 return D3(f, Nsteps, NE, step, i, j);
}