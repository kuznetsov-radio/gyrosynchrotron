import numpy as np
import matplotlib.pyplot as plt
import pygs  # Python wrapper located in the same directory. If not, need to include the path using sys.path.append()

# libname is the path to the C++ shared library. Modify as needed
libname = '../Binaries/MWTransferArr_OSX.so'

#
Nf = 100  # number of frequencies
NSteps = 10  # number of nodes along the line-of-sight
N_E = 15  # number of energy nodes
N_mu = 15  # number of pitch-angle nodes

Lparms = np.zeros(11, dtype='int')  # array of dimensions etc.
Lparms[0] = NSteps
Lparms[1] = Nf

Rparms = np.zeros(5)  # array of global floating-point parameters
Rparms[0] = 1e20  # area, cm^2
Rparms[1] = 1e9  # starting frequency to calculate spectrum, Hz
Rparms[2] = 0.02  # logarithmic step in frequency
Rparms[3] = 12  # f^C
Rparms[4] = 12  # f^WH

L = 1e10  # total source depth, cm

ParmLocal = np.zeros(24)  # array of voxel parameters - for a single voxel
ParmLocal[0] = L / NSteps  # voxel depth, cm
ParmLocal[1] = 3e7  # T_0, K
ParmLocal[2] = 3e9  # n_0 - thermal electron density, cm^{-3}
ParmLocal[3] = 180  # B - magnetic field, G
ParmLocal[6] = 3  # distribution over energy (PLW is chosen)
ParmLocal[7] = 1e6  # n_b - nonthermal electron density, cm^{-3}
ParmLocal[9] = 0.1  # E_min, MeV
ParmLocal[10] = 10.0  # E_max, MeV
ParmLocal[12] = 4.0  # \delta_1
ParmLocal[14] = 3  # distribution over pitch-angle (GLC is chosen)
ParmLocal[15] = 70  # loss-cone boundary, degrees
ParmLocal[16] = 0.2  # \Delta\mu

Parms = np.zeros((NSteps, 24))  # 2D array of input parameters - for multiple voxels
for i in range(NSteps):
    Parms[i] = ParmLocal  # most of the parameters are the same in all voxels
    Parms[i, 4] = 50.0 + 30.0 * i / (NSteps - 1)  # the viewing angle varies from 50 to 80 degrees along the LOS

# calculating the emission for analytical distribution (array -> on)
RL = pygs.get_mw(libname, Lparms, Rparms, Parms, 0., 0., 0.)

# plot the results
f = RL[:, 0]  # emission frequency, GHz
I_L = RL[:, 5]  # left-polarized emission intensity (as observed from the Earth), sfu
I_R = RL[:, 6]  # right-polarized emission intensity (as observed from the Earth), sfu
flux = I_L + I_R
pol = (I_L - I_R) / (I_L + I_R)

# Plot the results
plt.figure(figsize=(12, 5))
ax = plt.subplot(121)
# first plot, total intensity
ax.plot(f, flux)
ax.set_xlabel('Frequency (GHz)')
ax.set_ylabel('Flux (sfu)')
ax.set_title('Total Intensity')
ax.set_xscale('log')
ax.set_yscale('log')
ax = plt.subplot(122)
# second plot, circular polarization degree
ax.plot(f, pol)
ax.set_xlabel('Frequency (GHz)')
ax.set_ylabel('Polarization Degree')
ax.set_title('Circular Polarization Degree')
ax.set_xscale('log')
ax.set_yscale('linear')
ax.set_ylim([-1, 1])
plt.tight_layout()
plt.show()
