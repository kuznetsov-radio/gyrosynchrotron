pro Example_Array_SingleThread
 ;libname='MWTransferArr64.dll' ;name of the executable library
 ;libname='../Binaries/MWTransferArr_MacOS.so' ;name of the executable library
 libname='../MWTransferArr.so' ;name of the executable library

 Nf=100     ;number of frequencies
 NSteps=30L ;number of nodes along the line-of-sight
 
 N_E=15     ;number of energy nodes
 N_mu=15    ;number of pitch-angle nodes

 Lparms=lonarr(11) ;array of dimensions etc.
 Lparms[0]=Nsteps
 Lparms[1]=Nf
 Lparms[2]=N_E
 Lparms[3]=N_mu
 
 Rparms=dblarr(5) ;array of global floating-point parameters
 Rparms[0]=1e20 ;area, cm^2
 Rparms[1]=1e9  ;starting frequency to calculate spectrum, Hz
 Rparms[2]=0.02 ;logarithmic step in frequency
 Rparms[3]=12   ;f^C
 Rparms[4]=12   ;f^WH
 
 L=1e10 ;total source depth, cm
 
 ParmLocal=dblarr(24) ;array of voxel parameters - for a single voxel
 ParmLocal[0]=L/NSteps  ;voxel depth, cm
 ParmLocal[1]=3e7   ;T_0, K
 ParmLocal[2]=3e9   ;n_0 - thermal electron density, cm^{-3}
 ParmLocal[3]=180   ;B - magnetic field, G
 
 Parms=dblarr(24, NSteps) ;2D array of input parameters - for multiple voxels
 for i=0, NSteps-1 do begin
  Parms[*, i]=ParmLocal ;most of the parameters are the same in all voxels
  Parms[4, i]=50.0+30.0*i/(NSteps-1) ;the viewing angle varies from 50 to 80 degrees along the LOS
 endfor
 
 ;parameters of the electron distribution function 
 Emin=0.1d0 ;E_min, MeV
 Emax=10d0  ;E_max, MeV
 delta=4d0  ;delta_1
 n_b=1e6    ;n_b - nonthermal electron density, cm^{-3}
 mu_c=cos(!dpi*70/180) ;loss-cone boundary
 dmu_c=0.2d0 ;Delta_mu
 
 E_arr=exp(alog(Emin)+alog(Emax/Emin)*dindgen(N_E)/(N_E-1)) ;energy grid (logarithmically spaced)
 mu_arr=-1d0+2d0*dindgen(N_mu)/(N_mu-1) ;pitch-angle grid
 
 f0=dblarr(N_E, N_mu) ;2D distribution function array - for a single voxel
 
 ;computing the distribution function (equivalent to PLW & GLC)
 A=n_b/(2d0*!dpi)*(delta-1d0)/(Emin^(1d0-delta)-Emax^(1d0-delta))
 B=0.5d0/(mu_c+dmu_c*sqrt(!dpi)/2*erf((1d0-mu_c)/dmu_c))
 for i=0L, N_E-1 do for j=0L, N_mu-1 do begin
  amu=abs(mu_arr[j])
  f0[i, j]=A*B*E_arr[i]^(-delta)*((amu lt mu_c) ? 1d0 : exp(-((amu-mu_c)/dmu_c)^2))
 endfor
 
 f_arr=dblarr(N_E, N_mu, Nsteps) ;3D distribution function array - for multiple voxels
 for k=0L, Nsteps-1 do f_arr[*, *, k]=f0 ;electron distribution function is the same in all voxels
 
 ;calculating the emission for analytical distribution (array -> on)
 RL=dblarr(7, Nf)
 
 res=call_external(libname, 'GET_MW', Lparms, Rparms, Parms, E_arr, mu_arr, f_arr, RL)
 
 ;retrieving the results
 f=RL[0, *]
 I_L=RL[5, *]
 I_R=RL[6, *]

 ;--------------------------------------------
 ;plotting the results
 
 window, 3, title='Total intensity (array)'
 wset, 3
 plot, f, I_L+I_R, /xlog, /ylog, xtitle='Frequency, GHz', ytitle='Intensity, sfu'
     
 window, 4, title='Circular polarization degree (array)'
 wset, 4
 plot, f, (I_L-I_R)/(I_L+I_R), /xlog, xtitle='Frequency, GHz', ytitle='Polarization degree'
end
