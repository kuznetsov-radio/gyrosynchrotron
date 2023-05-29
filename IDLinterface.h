#pragma once

#define RpSize 6
#define InSize 24
#define OutSize 7

#define i_Npix 0
#define i_Nz 0
#define i_Nnu 1
#define i_NE 2
#define i_Nmu 3
#define i_Nnodes 4
#define i_matchKey 5
#define i_QoptKey 6
#define i_arrKeyG 7
#define i_logkey 8
#define i_PKkey 9
#define i_splinekey 10

#define i_S 0
#define i_nu0 1
#define i_dnu 2
#define i_nuCr 3
#define i_nuWH 4
#define i_dSun 5

#define i_dz 0
#define i_T0 1
#define i_n0 2
#define i_B 3
#define i_theta 4
#define i_EMflag 5
#define i_EId 6
#define i_nb 7
#define i_epskappa 8
#define i_Emin 9
#define i_Emax 10
#define i_Ebreak 11
#define i_delta1 12
#define i_delta2 13
#define i_muId 14
#define i_alphaC 15
#define i_dmu 16
#define i_a4 17
#define i_np 18
#define i_nH 19
#define i_nHe 20
#define i_arrKeyL 21
#define i_abcode 22

#define iRL_nu 0
#define iRL_Lw 1
#define iRL_Rw 2
#define iRL_Ls 3
#define iRL_Rs 4
#define iRL_Le 5
#define iRL_Re 6

#ifndef LINUX
extern "C" __declspec(dllexport) int GET_MW(int argc, void **argv);
extern "C" __declspec(dllexport) int GET_MW_SLICE(int argc, void **argv);
#else
extern "C" double GET_MW(int argc, void **argv);
extern "C" double GET_MW_SLICE(int argc, void **argv);
#endif