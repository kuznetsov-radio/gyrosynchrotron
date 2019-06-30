#include <Windows.h>
#include "Messages.h"
#include "ExtMath.h"
#include "IDLInterface.h"

int ProcessorNumber;

typedef struct
{
 int *Ndat;
 double *Parms;
 double *E_arr;
 double *mu_arr;
 double *f_arr;
 double *RL;

 HANDLE M;
 int *Done;
} PixList;

DWORD WINAPI SliceThread(LPVOID arg)
{
 #define OutSize 7

 PixList *PL=(PixList*)arg;

 int AllDone=0;

 while (!AllDone)
 {
  int Npix=PL->Ndat[0];
  int pix=-1;

  WaitForSingleObject(PL->M, INFINITE);
  for (int i=0; i<Npix; i++) if (!PL->Done[i])
  {
   pix=i;
   PL->Done[i]=1;
   break;
  }
  ReleaseMutex(PL->M);

  if (pix<0) AllDone=1;
  else
  {
   int Nsteps=PL->Ndat[1];
   int NE=PL->Ndat[2];
   int Nmu=PL->Ndat[3];
   int InSize=PL->Ndat[4];
   int Nnu=(int)D3(PL->Parms, InSize, Nsteps, 18, 0, 0);

   if (pix>0)
   {
	D3(PL->Parms, InSize, Nsteps, 15, 0, pix)=D3(PL->Parms, InSize, Nsteps, 15, 0, 0);
	D3(PL->Parms, InSize, Nsteps, 16, 0, pix)=D3(PL->Parms, InSize, Nsteps, 16, 0, 0);
	D3(PL->Parms, InSize, Nsteps, 18, 0, pix)=D3(PL->Parms, InSize, Nsteps, 18, 0, 0);
	for (int i=0; i<Nnu; i++) D3(PL->RL, OutSize, Nnu, 0, i, pix)=D3(PL->RL, OutSize, Nnu, 0, i, 0);
   }

   void *ARGV[6];
   ARGV[0]=(void*)(PL->Ndat+1);
   ARGV[1]=(void*)(PL->Parms+pix*Nsteps*InSize);
   ARGV[2]=(void*)(PL->E_arr);
   ARGV[3]=(void*)(PL->mu_arr);
   ARGV[4]=(void*)(PL->f_arr+pix*Nsteps*NE*Nmu);
   ARGV[5]=(void*)(PL->RL+pix*Nnu*OutSize);

   GET_MW(6, ARGV);
  }
 }

 return 0;
}

extern "C" __declspec(dllexport) double GET_MW_SLICE(int argc, void **argv)
{
 if (argc<6)
 {
  IDLmsg("GET_MW_SLICE error: not enough parameters in the function call.");
  return -1;
 }

 PixList PL;
 PL.Ndat=(int*)argv[0];
 PL.Parms=(double*)argv[1];
 PL.E_arr=(double*)argv[2];
 PL.mu_arr=(double*)argv[3];
 PL.f_arr=(double*)argv[4];
 PL.RL=(double*)argv[5];

 int Npix=PL.Ndat[0];
 int Nsteps=PL.Ndat[1];
 int InSize=PL.Ndat[4];

 if (Npix<=0)
 {
  IDLmsg("GET_MW_SLICE error: number of pixels must be positive.");
  return -2;
 }

 PL.Done=(int*)malloc(sizeof(int)*Npix);
 for (int i=0; i<Npix; i++) PL.Done[i]=0;

 int NThreads=(int)PL.D3(Parms, InSize, Nsteps, 24, 0, 0);
 if (NThreads<1) NThreads=(ProcessorNumber>1) ? ProcessorNumber-1 : 1;
 HANDLE *ThreadList=(HANDLE*)malloc(sizeof(HANDLE)*NThreads);
 
 PL.M=CreateMutexA(0, FALSE, 0);

 for (int i=0; i<NThreads; i++) ThreadList[i]=CreateThread(0, 0, SliceThread, &PL, 0, 0);
 WaitForMultipleObjects(NThreads, ThreadList, TRUE, INFINITE);

 for (int i=0; i<NThreads; i++) CloseHandle(ThreadList[i]);

 CloseHandle(PL.M);
  
 free(PL.Done);
 free(ThreadList);

 return 0;
}