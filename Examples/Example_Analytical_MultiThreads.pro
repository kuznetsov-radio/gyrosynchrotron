pro Example_Analytical_MultiThreads
 libname='MWTransferArr64.dll'
 ;libname='MWTransferArr.so'

 Nf=100     ;number of frequencies
 NSteps=30L ;number of nodes along the line-of-sight
 Npix=4L    ;number of pixels
 
 ParmLocal=dblarr(34) ;the array of the volume element parameters
 ParmLocal[0] =1e20  ;Area, cm^2
 ParmLocal[1] =1e10  ;Depth, cm (will be changed later!)
 ParmLocal[2] =3e7   ;T_0, K
 ParmLocal[3] =0.05  ;\eps (not used in this example)
 ParmLocal[4] =4.0   ;\kappa (not used in this example)
 ParmLocal[5] =16    ;number of integration nodes
 ParmLocal[6] =0.1   ;E_min, MeV
 ParmLocal[7] =10.0  ;E_max, MeV
 ParmLocal[8] =1.0   ;E_break, MeV (not used in this example)
 ParmLocal[9] =4.0   ;\delta_1
 ParmLocal[10]=6.0   ;\delta_2 (not used in this example)
 ParmLocal[11]=3e9   ;n_0 - thermal electron density, cm^{-3}
 ParmLocal[12]=1e6   ;n_b - nonthermal electron density, cm^{-3} (will be changed later)
 ParmLocal[13]=180   ;B - magnetic field, G
 ParmLocal[14]=60    ;theta - the viewing angle, degrees (will be changed later!)
 ParmLocal[15]=1e9   ;starting frequency to calculate spectrum, Hz
 ParmLocal[16]=0.02  ;logarithmic step in frequency
 ParmLocal[17]=3     ;distribution over energy (PLW is chosen)
 ParmLocal[18]=Nf    ;number of frequencies (specified above)
 ParmLocal[19]=3     ;distribution over pitch-angle (GLC is chosen)
 ParmLocal[20]=90    ;loss-cone boundary, degrees (will be changed later)
 ParmLocal[21]=0     ;beam direction (degrees) in GAU and SGA (not used in this example)
 ParmLocal[22]=0.2   ;\Delta\mu
 ParmLocal[23]=1     ;a_4 in SGA (not used in this example)
 ParmLocal[25]=12    ;f^C_cr
 ParmLocal[26]=12    ;f^WH_cr
 ParmLocal[27]=1     ;matching on
 ParmLocal[28]=2     ;Q-optimization on
 ;other parameters are zero by default => array is off
 
 Parms=dblarr(34, NSteps) ;the array of input parameters
 for i=0, NSteps-1 do begin
  Parms[*, i] =ParmLocal     
  Parms[1, i]/=NSteps          
  Parms[14, i]=50.0+30.0*i/(NSteps-1) ;The viewing angle varies from 50 to 80 degrees.
 endfor
 
 ParmsS=dblarr(34, Nsteps, Npix) ;the array of input parameters - for several pixels
 for i=0, Npix-1 do begin
  ParmsS[*, *, i]=Parms
  ParmsS[12, *, i]=1e6-3e5*i   ;electron density decreases from 1e6 to 1e5 cm^{-3}
  ParmsS[20, *, i]=90.0-25.0*i ;loss-cone angle decreases from 90 to 15 degrees
 endfor
 
 Ndat=lonarr(5)
 Ndat[0]=Npix
 Ndat[1]=NSteps
 Ndat[2]=0
 Ndat[3]=0
 Ndat[4]=34
 
 ;--------------------------------------------
 
 RL=dblarr(7, Nf, Npix)
 
 res=call_external(libname, 'GET_MW_SLICE', Ndat, ParmsS, 0, 0, 0, RL, /d_value)

 ;--------------------------------------------
 
 window, 1, title='Total intensity (analytical)'
 wset, 1
 
 for i=0, Npix-1 do begin
  f=RL[0, *, i]
  I_L=RL[5, *, i]
  I_R=RL[6, *, i]
  if i eq 0 then plot, f, I_L+I_R, /xlog, /ylog, $
                       xtitle='Frequency, GHz', ytitle='Intensity, sfu' $
                 else oplot, f, I_L+I_R, linestyle=i
 endfor                
   
 window, 2, title='Circular polarization degree (analytical)'
 wset, 2
 
 for i=0, Npix-1 do begin
  f=RL[0, *, i]
  I_L=RL[5, *, i]
  I_R=RL[6, *, i]
  if i eq 0 then plot, f, (I_L-I_R)/(I_L+I_R), /xlog, $
                       xtitle='Frequency, GHz', ytitle='Polarization degree', $
                       yrange=[-0.7, 0.7], ystyle=1 $
            else oplot, f, (I_L-I_R)/(I_L+I_R), linestyle=i
 endfor
end
