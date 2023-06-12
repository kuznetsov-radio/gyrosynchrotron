import numpy as np
import math
import matplotlib.pyplot as plt
import GScodes # initialization library - located either in the current directory or in the system path

libname='./MWTransferArr64.dll' # name of the executable library - located where Python can find it

GET_MW_SLICE=GScodes.initGET_MW_SLICE(libname) # load the library

Nf=100     # number of frequencies
NSteps=30  # number of nodes along the line-of-sight
Npix=4     # number of pixels
 
N_E=15     # number of energy nodes
N_mu=15    # number of pitch-angle nodes
 
Lparms_M=np.zeros(12, dtype='int32') # array of dimensions etc.
Lparms_M[0]=Npix
Lparms_M[1]=NSteps
Lparms_M[2]=Nf
Lparms_M[3]=N_E
Lparms_M[4]=N_mu
 
Rparms=np.zeros(5, dtype='double') # array of global floating-point parameters - for a single LOS
Rparms[0]=1e20 # area, cm^2
Rparms[1]=1e9  # starting frequency to calculate spectrum, Hz
Rparms[2]=0.02 # logarithmic step in frequency
Rparms[3]=12   # f^C_cr
Rparms[4]=12   # f^WH_cr
 
Rparms_M=np.zeros((5, Npix), dtype='double', order='F') # 2D array of global floating-point parameters - for multiple LOSs
for pix in range(Npix):
    Rparms_M[:, pix]=Rparms # the parameters are the same for all LOSs
 
L=1e10 # source depth, cm
 
ParmLocal=np.zeros(24, dtype='double') # array of voxel parameters - for a single voxel
ParmLocal[0]=L/NSteps  # voxel depth, cm
ParmLocal[1]=3e7   # T_0, K
ParmLocal[2]=3e9   # n_0 - thermal electron density, cm^{-3}
ParmLocal[3]=180   # B - magnetic field, G
 
Parms=np.zeros((24, NSteps), dtype='double', order='F') # 2D array of input parameters - for multiple voxels
for i in range(NSteps):
    Parms[:, i]=ParmLocal # most of the parameters are the same in all voxels
    Parms[4, i]=50.0+30.0*i/(NSteps-1) # the viewing angle varies from 50 to 80 degrees along the LOS
 
Parms_M=np.zeros((24, NSteps, Npix), dtype='double', order='F') # 3D array of input parameters - for multiple voxels and LOSs
for pix in range(Npix):
    Parms_M[:, :, pix]=Parms # parameters are the same for all LOSs
 
# common parameters of the electron distribution function 
Emin=0.1  # E_min, MeV
Emax=10.0 # E_max, MeV
delta=4.0 # delta_1
dmu_c=0.2 # Delta_mu
 
E_arr=np.logspace(np.log10(Emin), np.log10(Emax), N_E, dtype='double') # energy grid (logarithmically spaced)
mu_arr=np.linspace(-1.0, 1.0, N_mu, dtype='double') # pitch-angle grid
 
f_arr_M=np.zeros((N_E, N_mu, NSteps, Npix), dtype='double', order='F') # 4D distribution function array
 
for k in range(Npix):
    n_b=1e6-3e5*k  # electron density decreases from 1e6 to 1e5 cm^{-3}, depending on the LOS
    mu_c=np.cos((90.0-25.0*k)*np.pi/180) # loss-cone angle decreases from 90 to 15 degrees, depending on the LOS
  
    f0=np.zeros((N_E, N_mu), dtype='double') # 2D distribution function array - for a single voxel
 
    # computing the distribution function (equivalent to PLW & GLC)
    A=n_b/(2.0*np.pi)*(delta-1.0)/(Emin**(1.0-delta)-Emax**(1.0-delta))
    B=0.5/(mu_c+dmu_c*np.sqrt(np.pi)/2*math.erf((1.0-mu_c)/dmu_c))
    for i in range(N_E):
        for j in range(N_mu):
            amu=abs(mu_arr[j])
            f0[i, j]=A*B*E_arr[i]**(-delta)*(1.0 if amu<mu_c else np.exp(-((amu-mu_c)/dmu_c)**2))

    for l in range(NSteps):
        f_arr_M[:, :, l, k]=f0 # electron distribution function is the same for all voxels within a LOS, but varies from LOS to LOS
 
RL_M=np.zeros((7, Nf, Npix), dtype='double', order='F') # input/output array

# calculating the emission for array distribution (array -> on) 
res=GET_MW_SLICE(Lparms_M, Rparms_M, Parms_M, E_arr, mu_arr, f_arr_M, RL_M)
 
# plotting the results
 
plt.figure(1)
for i in range(Npix):
    f=RL_M[0, :, i]
    I_L=RL_M[5, :, i]
    I_R=RL_M[6, :, i]
    plt.plot(f, I_L+I_R)
plt.xscale('log')
plt.yscale('log')
plt.title('Total intensity (array)')
plt.xlabel('Frequency, GHz')
plt.ylabel('Intensity, sfu')       
   
plt.figure(2)
for i in range(Npix):
    f=RL_M[0, :, i]
    I_L=RL_M[5, :, i]
    I_R=RL_M[6, :, i]
    plt.plot(f, (I_L-I_R)/(I_L+I_R))
plt.xscale('log')
plt.title('Circular polarization degree (array)')
plt.xlabel('Frequency, GHz')
plt.ylabel('Polarization degree')

plt.show()    
