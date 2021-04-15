pro Example_Array_MultiThreads
 libname='MWTransferArr64.dll' ;name of the executable library

 Nf=100     ;number of frequencies
 NSteps=30L ;number of nodes along the line-of-sight
 Npix=4L    ;number of pixels
 
 N_E=15     ;number of energy nodes
 N_mu=15    ;number of pitch-angle nodes
 
 Lparms_M=lonarr(12) ;array of dimensions etc.
 Lparms_M[0]=Npix
 Lparms_M[1]=Nsteps
 Lparms_M[2]=Nf
 Lparms_M[3]=N_E
 Lparms_M[4]=N_mu
 
 Rparms=dblarr(5) ;array of global floating-point parameters - for a single LOS
 Rparms[0]=1e20 ;area, cm^2
 Rparms[1]=1e9  ;starting frequency to calculate spectrum, Hz
 Rparms[2]=0.02 ;logarithmic step in frequency
 Rparms[3]=12   ;f^C_cr
 Rparms[4]=12   ;f^WH_cr
 
 Rparms_M=dblarr(5, Npix) ;2D array of global floating-point parameters - for multiple LOSs
 for pix=0, Npix-1 do Rparms_M[*, pix]=Rparms ;the parameters are the same for all LOSs
 
 L=1e10 ;source depth, cm
 
 ParmLocal=dblarr(24) ;array of voxel parameters - for a single voxel
 ParmLocal[0]=L/NSteps  ;voxel depth, cm
 ParmLocal[1]=3e7   ;T_0, K
 ParmLocal[2]=3e9   ;n_0 - thermal electron density, cm^{-3}
 ParmLocal[3]=180   ;B - magnetic field, G
 
 Parms=dblarr(24, NSteps) ;2D array of input parameters - for multiple voxels
 for i=0, NSteps-1 do begin
  Parms[*, i] =ParmLocal ;most of the parameters are the same in all voxels
  Parms[4, i]=50.0+30.0*i/(NSteps-1) ;the viewing angle varies from 50 to 80 degrees along the LOS
 endfor
 
 Parms_M=dblarr(24, Nsteps, Npix) ;3D array of input parameters - for multiple voxels and LOSs
 for pix=0, Npix-1 do Parms_M[*, *, pix]=Parms ;parameters are the same for all LOSs
 
 ;common parameters of the electron distribution function 
 Emin=0.1d0 ;E_min, MeV
 Emax=10d0  ;E_max, MeV
 delta=4d0  ;delta_1
 dmu_c=0.2d0 ;Delta_mu
 
 E_arr=exp(alog(Emin)+alog(Emax/Emin)*dindgen(N_E)/(N_E-1)) ;energy grid (logarithmically spaced)
 mu_arr=-1d0+2d0*dindgen(N_mu)/(N_mu-1) ;pitch-angle grid
 
 f_arr_M=dblarr(N_E, N_mu, NSteps, Npix) ;4D distribution function array
 
 for k=0, Npix-1 do begin
  n_b=1d6-3d5*k  ;electron density decreases from 1e6 to 1e5 cm^{-3}, depending on the LOS
  mu_c=cos((90.0-25.0*k)*!dpi/180) ;loss-cone angle decreases from 90 to 15 degrees, depending on the LOS
  
  f0=dblarr(N_E, N_mu) ;2D distribution function array - for a single voxel
 
  ;computing the distribution function (equivalent to PLW & GLC)
  A=n_b/(2d0*!dpi)*(delta-1d0)/(Emin^(1d0-delta)-Emax^(1d0-delta))
  B=0.5d0/(mu_c+dmu_c*sqrt(!dpi)/2*erf((1d0-mu_c)/dmu_c))
  for i=0L, N_E-1 do for j=0L, N_mu-1 do begin
   amu=abs(mu_arr[j])
   f0[i, j]=A*B*E_arr[i]^(-delta)*((amu lt mu_c) ? 1d0 : exp(-((amu-mu_c)/dmu_c)^2))
  endfor
   
  for l=0, NSteps-1 do f_arr_M[*, *, l, k]=f0 
  ;electron distribution function is the same for all voxels within a LOS, but varies from LOS to LOS
 endfor 
 
 ;calculating the emission for analytical distribution (array -> on)
 RL_M=dblarr(7, Nf, Npix)
 
 res=call_external(libname, 'GET_MW_SLICE', Lparms_M, Rparms_M, Parms_M, E_arr, mu_arr, f_arr_M, RL_M, /unload)
 
 ;--------------------------------------------
 ;plotting the results
 
 window, 3, title='Total intensity (array)'
 wset, 3
  
 for i=0, Npix-1 do begin
  f=RL_M[0, *, i]
  I_L=RL_M[5, *, i]
  I_R=RL_M[6, *, i]
  if i eq 0 then plot, f, I_L+I_R, /xlog, /ylog, $
                       xtitle='Frequency, GHz', ytitle='Intensity, sfu' $
                 else oplot, f, I_L+I_R, linestyle=i
 endfor                
   
 window, 4, title='Circular polarization degree (array)'
 wset, 4
 
 for i=0, Npix-1 do begin
  f=RL_M[0, *, i]
  I_L=RL_M[5, *, i]
  I_R=RL_M[6, *, i]
  if i eq 0 then plot, f, (I_L-I_R)/(I_L+I_R), /xlog, $
                       xtitle='Frequency, GHz', ytitle='Polarization degree', $
                       yrange=[-0.7, 0.7], ystyle=1 $
            else oplot, f, (I_L-I_R)/(I_L+I_R), linestyle=i
 endfor
end
