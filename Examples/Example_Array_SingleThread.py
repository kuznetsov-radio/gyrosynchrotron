import numpy as np
import math
import matplotlib.pyplot as plt
import GScodes # initialization library - located either in the current directory or in the system path

libname='./MWTransferArr64.dll' # name of the executable library - located where Python can find it

GET_MW=GScodes.initGET_MW(libname) # load the library

Nf=100     # number of frequencies
NSteps=30  # number of nodes along the line-of-sight
 
N_E=15     # number of energy nodes
N_mu=15    # number of pitch-angle nodes

Lparms=np.zeros(11, dtype='int32') # array of dimensions etc.
Lparms[0]=NSteps
Lparms[1]=Nf
Lparms[2]=N_E
Lparms[3]=N_mu
 
Rparms=np.zeros(5, dtype='double') # array of global floating-point parameters
Rparms[0]=1e20 # area, cm^2
Rparms[1]=1e9  # starting frequency to calculate spectrum, Hz
Rparms[2]=0.02 # logarithmic step in frequency
Rparms[3]=12   # f^C
Rparms[4]=12   # f^WH
 
L=1e10 # total source depth, cm
 
ParmLocal=np.zeros(24, dtype='double') # array of voxel parameters - for a single voxel
ParmLocal[0]=L/NSteps  # voxel depth, cm
ParmLocal[1]=3e7   # T_0, K
ParmLocal[2]=3e9   # n_0 - thermal electron density, cm^{-3}
ParmLocal[3]=180   # B - magnetic field, G
 
Parms=np.zeros((24, NSteps), dtype='double', order='F') # 2D array of input parameters - for multiple voxels
for i in range(NSteps-1):
    Parms[:, i]=ParmLocal # most of the parameters are the same in all voxels
    Parms[4, i]=50.0+30.0*i/(NSteps-1) # the viewing angle varies from 50 to 80 degrees along the LOS
 
# parameters of the electron distribution function 
Emin=0.1   # E_min, MeV
Emax=10.0  # E_max, MeV
delta=4.0  # delta_1
n_b=1e6    # n_b - nonthermal electron density, cm^{-3}
mu_c=np.cos(np.pi*70/180) # loss-cone boundary
dmu_c=0.2  # Delta_mu
 
E_arr=np.logspace(np.log10(Emin), np.log10(Emax), N_E, dtype='double') # energy grid (logarithmically spaced)
mu_arr=np.linspace(-1.0, 1.0, N_mu, dtype='double') # pitch-angle grid
 
f0=np.zeros((N_E, N_mu), dtype='double') # 2D distribution function array - for a single voxel
 
# computing the distribution function (equivalent to PLW & GLC)
A=n_b/(2.0*np.pi)*(delta-1.0)/(Emin**(1.0-delta)-Emax**(1.0-delta))
B=0.5/(mu_c+dmu_c*np.sqrt(np.pi)/2*math.erf((1.0-mu_c)/dmu_c))
for i in range(N_E):
    for j in range(N_mu):
        amu=abs(mu_arr[j])
        f0[i, j]=A*B*E_arr[i]**(-delta)*(1.0 if amu<mu_c else np.exp(-((amu-mu_c)/dmu_c)**2))
 
f_arr=np.zeros((N_E, N_mu, NSteps), dtype='double', order='F') # 3D distribution function array - for multiple voxels
for k in range(NSteps):
    f_arr[:, :, k]=f0 # electron distribution function is the same in all voxels
 
RL=np.zeros((7, Nf), dtype='double', order='F') # input/output array

# calculating the emission for array distribution (array -> on) 
res=GET_MW(Lparms, Rparms, Parms, E_arr, mu_arr, f_arr, RL)
 
# retrieving the results
f=RL[0]
I_L=RL[5]
I_R=RL[6]

# plotting the results
plt.figure(1)
plt.plot(f, I_L+I_R)
plt.xscale('log')
plt.yscale('log')
plt.title('Total intensity (array)')
plt.xlabel('Frequency, GHz')
plt.ylabel('Intensity, sfu')

plt.figure(2)
plt.plot(f, (I_L-I_R)/(I_L+I_R))
plt.xscale('log')
plt.title('Circular polarization degree (array)')
plt.xlabel('Frequency, GHz')
plt.ylabel('Polarization degree')

plt.show()
