pro Example_Analytical_SingleThread
 libname='MWTransferArr64.dll'

 Nf=100     ;number of frequencies
 NSteps=30L ;number of nodes along the line-of-sight
 
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
 ParmLocal[12]=1e6   ;n_b - nonthermal electron density, cm^{-3}
 ParmLocal[13]=180   ;B - magnetic field, G
 ParmLocal[14]=60    ;theta - the viewing angle, degrees (will be changed later!)
 ParmLocal[15]=1e9   ;starting frequency to calculate spectrum, Hz
 ParmLocal[16]=0.02  ;logarithmic step in frequency
 ParmLocal[17]=3     ;distribution over energy (PLW is chosen)
 ParmLocal[18]=Nf    ;number of frequencies (specified above)
 ParmLocal[19]=3     ;distribution over pitch-angle (GAU is chosen)
 ParmLocal[20]=70    ;loss-cone boundary, degrees
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
 
 Ndat=lonarr(4)
 Ndat[0]=NSteps
 Ndat[1]=0
 Ndat[2]=0
 Ndat[3]=34
 
 ;--------------------------------------------
 
 ;calculating the emission for analytical distribution (array -> off)
 RL=dblarr(7, Nf)
 
 res=call_external(libname, 'GET_MW', Ndat, Parms, 0, 0, 0, RL, /d_value)
 
 f=RL[0, *]
 I_L=RL[5, *]
 I_R=RL[6, *]
 
 ;--------------------------------------------
 
 set_plot, 'win'
  
 window, 1, title='Total intensity (analytical)'
 wset, 1
 plot, f, I_L+I_R, /xlog, /ylog, xtitle='Frequency, GHz', ytitle='Intensity, sfu'
    
 window, 2, title='Circular polarization degree (analytical)'
 wset, 2
 plot, f, (I_L-I_R)/(I_L+I_R), /xlog, xtitle='Frequency, GHz', ytitle='Polarization degree'
end