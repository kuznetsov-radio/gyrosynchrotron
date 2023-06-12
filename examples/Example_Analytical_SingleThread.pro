pro Example_Analytical_SingleThread
 libname='MWTransferArr64.dll' ;name of the executable library

 Nf=100     ;number of frequencies
 NSteps=30L ;number of nodes along the line-of-sight
 
 Lparms=lonarr(11) ;array of dimensions etc.
 Lparms[0]=Nsteps
 Lparms[1]=Nf
 
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
 ParmLocal[6]=3     ;distribution over energy (PLW is chosen)
 ParmLocal[7]=1e6   ;n_b - nonthermal electron density, cm^{-3}
 ParmLocal[9]=0.1   ;E_min, MeV
 ParmLocal[10]=10.0 ;E_max, MeV
 ParmLocal[12]=4.0  ;\delta_1
 ParmLocal[14]=3    ;distribution over pitch-angle (GLC is chosen)
 ParmLocal[15]=70   ;loss-cone boundary, degrees
 ParmLocal[16]=0.2  ;\Delta\mu
 
 Parms=dblarr(24, NSteps) ;2D array of input parameters - for multiple voxels
 for i=0, NSteps-1 do begin
  Parms[*, i]=ParmLocal ;most of the parameters are the same in all voxels
  Parms[4, i]=50.0+30.0*i/(NSteps-1) ;the viewing angle varies from 50 to 80 degrees along the LOS
 endfor
 
 ;calculating the emission for analytical distribution (array -> off)
 RL=dblarr(7, Nf)
 
 res=call_external(libname, 'GET_MW', Lparms, Rparms, Parms, 0, 0, 0, RL)
 
 ;retrieving the results
 f=RL[0, *]
 I_L=RL[5, *]
 I_R=RL[6, *]
 
 ;--------------------------------------------
 ;plotting the results
 
 window, 1, title='Total intensity (analytical)'
 wset, 1
 plot, f, I_L+I_R, /xlog, /ylog, xtitle='Frequency, GHz', ytitle='Intensity, sfu'
    
 window, 2, title='Circular polarization degree (analytical)'
 wset, 2
 plot, f, (I_L-I_R)/(I_L+I_R), /xlog, xtitle='Frequency, GHz', ytitle='Polarization degree'
end
