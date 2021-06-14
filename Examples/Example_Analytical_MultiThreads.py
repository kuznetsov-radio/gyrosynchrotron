import numpy as np
import matplotlib.pyplot as plt
import GScodes # initialization library - located either in the current directory or in the system path

libname='./MWTransferArr64.dll' # name of the executable library - located where Python can find it

GET_MW_SLICE=GScodes.initGET_MW_SLICE(libname) # load the library

Nf=100     # number of frequencies
NSteps=30  # number of nodes along the line-of-sight
Npix=4     # number of pixels
 
Lparms_M=np.zeros(12, dtype='int32') # array of dimensions etc.
Lparms_M[0]=Npix
Lparms_M[1]=NSteps
Lparms_M[2]=Nf
 
Rparms=np.zeros(5, dtype='double') # array of global floating-point parameters - for a single LOS
Rparms[0]=1e20 # area, cm^2
Rparms[1]=1e9  # starting frequency to calculate spectrum, Hz
Rparms[2]=0.02 # logarithmic step in frequency
Rparms[3]=12   # f^C
Rparms[4]=12   # f^WH
 
Rparms_M=np.zeros((5, Npix), dtype='double', order='F') # 2D array of global floating-point parameters - for multiple LOSs
for pix in range(Npix):
    Rparms_M[:, pix]=Rparms # the parameters are the same for all LOSs
 
L=1e10 # total source depth, cm
 
ParmLocal=np.zeros(24, dtype='double') # array of voxel parameters - for a single voxel
ParmLocal[0]=L/NSteps  # voxel depth, cm
ParmLocal[1]=3e7   # T_0, K
ParmLocal[2]=3e9   # n_0 - thermal electron density, cm^{-3}
ParmLocal[3]=180   # B - magnetic field, G
ParmLocal[6]=3     # distribution over energy (PLW is chosen)
ParmLocal[9]=0.1   # E_min, MeV
ParmLocal[10]=10.0 # E_max, MeV
ParmLocal[12]=4.0  # \delta_1
ParmLocal[14]=3    # distribution over pitch-angle (GLC is chosen)
ParmLocal[16]=0.2  # \Delta\mu
 
Parms=np.zeros((24, NSteps), dtype='double', order='F') # 2D array of input parameters - for multiple voxels
for i in range(NSteps):
    Parms[:, i]=ParmLocal # most of the parameters are the same in all voxels     
    Parms[4, i]=50.0+30.0*i/(NSteps-1) # the viewing angle varies from 50 to 80 degrees along the LOS
 
Parms_M=np.zeros((24, NSteps, Npix), dtype='double', order='F') # 3D array of input parameters - for multiple voxels and LOSs
for pix in range(Npix):
    Parms_M[:, :, pix]=Parms # most of the parameters are the same for all LOSs
    Parms_M[7, :, pix]=1e6-3e5*pix # electron density decreases from 1e6 to 1e5 cm^{-3}, depending on the LOS
    Parms_M[15, :, pix]=90.0-25.0*pix # loss-cone angle decreases from 90 to 15 degrees, depending on the LOS

RL_M=np.zeros((7, Nf, Npix), dtype='double', order='F') # input/output array
dummy=np.array(0, dtype='double')

# calculating the emission for analytical distribution (array -> off),
# the unused parameters can be set to any value
res=GET_MW_SLICE(Lparms_M, Rparms_M, Parms_M, dummy, dummy, dummy, RL_M)
 
# plotting the results
 
plt.figure(1)
for i in range(Npix):
    f=RL_M[0, :, i]
    I_L=RL_M[5, :, i]
    I_R=RL_M[6, :, i]
    plt.plot(f, I_L+I_R)
plt.xscale('log')
plt.yscale('log')
plt.title('Total intensity (analytical)')
plt.xlabel('Frequency, GHz')
plt.ylabel('Intensity, sfu')       
   
plt.figure(2)
for i in range(Npix):
    f=RL_M[0, :, i]
    I_L=RL_M[5, :, i]
    I_R=RL_M[6, :, i]
    plt.plot(f, (I_L-I_R)/(I_L+I_R))
plt.xscale('log')
plt.title('Circular polarization degree (analytical)')
plt.xlabel('Frequency, GHz')
plt.ylabel('Polarization degree')

plt.show()    
