#include "Messages.h"
#include "MWmain.h"
#include "ExtMath.h"
#include "IDLInterface.h"

extern "C" __declspec(dllexport) double GET_MW(int argc, void **argv)
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