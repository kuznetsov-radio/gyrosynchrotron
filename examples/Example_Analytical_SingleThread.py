import numpy as np
import matplotlib.pyplot as plt
import GScodes # initialization library - located either in the current directory or in the system path

libname='./MWTransferArr64.dll' # name of the executable library - located where Python can find it

GET_MW=GScodes.initGET_MW(libname) # load the library

Nf=100     # number of frequencies
NSteps=30  # number of nodes along the line-of-sight
 
Lparms=np.zeros(11, dtype='int32') # array of dimensions etc.
Lparms[0]=NSteps
Lparms[1]=Nf
 
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
ParmLocal[6]=3     # distribution over energy (PLW is chosen)
ParmLocal[7]=1e6   # n_b - nonthermal electron density, cm^{-3}
ParmLocal[9]=0.1   # E_min, MeV
ParmLocal[10]=10.0 # E_max, MeV
ParmLocal[12]=4.0  # \delta_1
ParmLocal[14]=3    # distribution over pitch-angle (GLC is chosen)
ParmLocal[15]=70   # loss-cone boundary, degrees
ParmLocal[16]=0.2  # \Delta\mu
 
Parms=np.zeros((24, NSteps), dtype='double', order='F') # 2D array of input parameters - for multiple voxels
for i in range(NSteps):
    Parms[:, i]=ParmLocal # most of the parameters are the same in all voxels
    Parms[4, i]=50.0+30.0*i/(NSteps-1) # the viewing angle varies from 50 to 80 degrees along the LOS
 
RL=np.zeros((7, Nf), dtype='double', order='F') # input/output array
dummy=np.array(0, dtype='double')

# calculating the emission for analytical distribution (array -> off),
# the unused parameters can be set to any value
res=GET_MW(Lparms, Rparms, Parms, dummy, dummy, dummy, RL)
 
# retrieving the results
f=RL[0]
I_L=RL[5]
I_R=RL[6]

# plotting the results
plt.figure(1)
plt.plot(f, I_L+I_R)
plt.xscale('log')
plt.yscale('log')
plt.title('Total intensity (analytical)')
plt.xlabel('Frequency, GHz')
plt.ylabel('Intensity, sfu')

plt.figure(2)
plt.plot(f, (I_L-I_R)/(I_L+I_R))
plt.xscale('log')
plt.title('Circular polarization degree (analytical)')
plt.xlabel('Frequency, GHz')
plt.ylabel('Polarization degree')

plt.show()
