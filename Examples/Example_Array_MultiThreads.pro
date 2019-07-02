pro Example_Array_MultiThreads
 libname='MWTransferArr64.dll'
 ;libname='MWTransferArr.so'

 Nf=100     ;number of frequencies
 NSteps=30L ;number of nodes along the line-of-sight
 Npix=4L    ;number of pixels
 
 N_E=15     ;number of energy nodes
 N_mu=15    ;number of pitch-angle nodes
 
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
 ParmLocal[17]=1     ;distribution over energy (free-free only)
 ParmLocal[18]=Nf    ;number of frequencies (specified above)
 ParmLocal[19]=3     ;distribution over pitch-angle (GLC is chosen)
 ParmLocal[20]=90    ;loss-cone boundary, degrees (will be changed later)
 ParmLocal[21]=0     ;beam direction (degrees) in GAU and SGA (not used in this example)
 ParmLocal[22]=0.2   ;\Delta\mu
 ParmLocal[23]=1     ;a_4 in SGA (not used in this example)
 ParmLocal[25]=12     ;f^C_cr
 ParmLocal[26]=12     ;f^WH_cr
 ParmLocal[27]=1     ;matching on
 ParmLocal[28]=2     ;Q-optimization on
 ParmLocal[31]=1     ;array on
 ParmLocal[32]=2     ;full array on 
 
 Parms=dblarr(34, NSteps) ;the array of input parameters
 for i=0, NSteps-1 do begin
  Parms[*, i] =ParmLocal     
  Parms[1, i]/=NSteps          
  Parms[14, i]=50.0+30.0*i/(NSteps-1) ;The viewing angle varies from 50 to 80 degrees.
 endfor
 
 ParmsS=dblarr(34, Nsteps, Npix) ;the array of input parameters - for several pixels
 for i=0, Npix-1 do ParmsS[*, *, i]=Parms
 
 Ndat=lonarr(5)
 Ndat[0]=Npix
 Ndat[1]=NSteps
 Ndat[2]=N_E
 Ndat[3]=N_mu
 Ndat[4]=34
 
 Emin=0.1d0 ;E_min, MeV
 Emax=10d0  ;E_max, MeV
 delta=4d0  ;delta_1
 dmu_c=0.2d0 ;Delta_mu
 
 E=exp(alog(Emin)+alog(Emax/Emin)*dindgen(N_E)/(N_E-1))
 mu=-1d0+2d0*dindgen(N_mu)/(N_mu-1)
 
 f=dblarr(Nsteps, N_E, N_mu, Npix)
 
 for k=0, Npix-1 do begin
  n_b=1d6-3d5*k   ;electron density decreases from 1e6 to 1e5 cm^{-3}
  mu_c=cos((90.0-25.0*k)*!dpi/180) ;loss-cone angle decreases from 90 to 15 degrees
  
  f0=dblarr(N_E, N_mu)
 
  A=n_b/(2d0*!dpi)*(delta-1d0)/(Emin^(1d0-delta)-Emax^(1d0-delta))
  B=0.5d0/(mu_c+dmu_c*sqrt(!dpi)/2*erf((1d0-mu_c)/dmu_c))
  for i=0L, N_E-1 do for j=0L, N_mu-1 do begin
   amu=abs(mu[j])
   f0[i, j]=A*B*E[i]^(-delta)*((amu lt mu_c) ? 1d0 : exp(-((amu-mu_c)/dmu_c)^2))
  endfor
   
  for i=0L, N_E-1 do for j=0L, N_mu-1 do f[*, i, j, k]=f0[i, j]
 endfor 
 
 ;--------------------------------------------
 
 RL=dblarr(7, Nf, Npix)
 
 res=call_external(libname, 'GET_MW_SLICE', Ndat, ParmsS, E, mu, f, RL, /d_value)
 
 ;--------------------------------------------
 
 window, 3, title='Total intensity (array)'
 wset, 3
 
 for i=0, Npix-1 do begin
  f=RL[0, *, i]
  I_L=RL[5, *, i]
  I_R=RL[6, *, i]
  if i eq 0 then plot, f, I_L+I_R, /xlog, /ylog, $
                       xtitle='Frequency, GHz', ytitle='Intensity, sfu' $
                 else oplot, f, I_L+I_R, linestyle=i
 endfor                
   
 window, 4, title='Circular polarization degree (array)'
 wset, 4
 
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
