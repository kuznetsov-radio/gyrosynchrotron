pro Example_Analytical_MultiThreads
 libname='MWTransferArr64.dll' ;name of the executable library

 Nf=100     ;number of frequencies
 NSteps=30L ;number of nodes along the line-of-sight
 Npix=4L    ;number of pixels
 
 Lparms_M=lonarr(12) ;array of dimensions etc.
 Lparms_M[0]=Npix
 Lparms_M[1]=Nsteps
 Lparms_M[2]=Nf
 
 Rparms=dblarr(5) ;array of global floating-point parameters - for a single LOS
 Rparms[0]=1e20 ;area, cm^2
 Rparms[1]=1e9  ;starting frequency to calculate spectrum, Hz
 Rparms[2]=0.02 ;logarithmic step in frequency
 Rparms[3]=12   ;f^C
 Rparms[4]=12   ;f^WH
 
 Rparms_M=dblarr(5, Npix) ;2D array of global floating-point parameters - for multiple LOSs
 for pix=0, Npix-1 do Rparms_M[*, pix]=Rparms ;the parameters are the same for all LOSs
 
 L=1e10 ;total source depth, cm
 
 ParmLocal=dblarr(24) ;array of voxel parameters - for a single voxel
 ParmLocal[0]=L/NSteps  ;voxel depth, cm
 ParmLocal[1]=3e7   ;T_0, K
 ParmLocal[2]=3e9   ;n_0 - thermal electron density, cm^{-3}
 ParmLocal[3]=180   ;B - magnetic field, G
 ParmLocal[6]=3     ;distribution over energy (PLW is chosen)
 ParmLocal[9]=0.1   ;E_min, MeV
 ParmLocal[10]=10.0 ;E_max, MeV
 ParmLocal[12]=4.0  ;\delta_1
 ParmLocal[14]=3    ;distribution over pitch-angle (GLC is chosen)
 ParmLocal[16]=0.2  ;\Delta\mu
 
 Parms=dblarr(24, NSteps) ;2D array of input parameters - for multiple voxels
 for i=0, NSteps-1 do begin
  Parms[*, i]=ParmLocal ;most of the parameters are the same in all voxels     
  Parms[4, i]=50.0+30.0*i/(NSteps-1) ;the viewing angle varies from 50 to 80 degrees along the LOS
 endfor
 
 Parms_M=dblarr(24, Nsteps, Npix) ;3D array of input parameters - for multiple voxels and LOSs
 for pix=0, Npix-1 do begin
  Parms_M[*, *, pix]=Parms ;most of the parameters are the same for all LOSs
  Parms_M[7, *, pix]=1e6-3e5*pix ;electron density decreases from 1e6 to 1e5 cm^{-3}, depending on the LOS
  Parms_M[15, *, pix]=90.0-25.0*pix ;loss-cone angle decreases from 90 to 15 degrees, depending on the LOS
 endfor

 ;calculating the emission for analytical distribution (array -> off)
 RL_M=dblarr(7, Nf, Npix)
 
 res=call_external(libname, 'GET_MW_SLICE', Lparms_M, Rparms_M, Parms_M, 0, 0, 0, RL_M, /unload)
 
 ;--------------------------------------------
 ;plotting the results
 
 window, 1, title='Total intensity (analytical)'
 wset, 1
 
 for i=0, Npix-1 do begin
  f=RL_M[0, *, i]
  I_L=RL_M[5, *, i]
  I_R=RL_M[6, *, i]
  if i eq 0 then plot, f, I_L+I_R, /xlog, /ylog, $
                       xtitle='Frequency, GHz', ytitle='Intensity, sfu' $
                 else oplot, f, I_L+I_R, linestyle=i
 endfor                
   
 window, 2, title='Circular polarization degree (analytical)'
 wset, 2
 
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
