import numpy as np
import math
import matplotlib.pyplot as plt
import pygs  # Python wrapper located in the same directory. If not, need to include the path using sys.path.append()

# libname is the path to the C++ shared library. Modify as needed
libname = '../Binaries/MWTransferArr.so'

#
Nf = 100  # number of frequencies
NSteps = 10  # number of nodes along the line-of-sight
N_E = 15  # number of energy nodes
N_mu = 15  # number of pitch-angle nodes

Lparms = np.zeros(11, dtype='int')  # array of dimensions etc.
Lparms[0] = NSteps
Lparms[1] = Nf
Lparms[2] = N_E
Lparms[3] = N_mu

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

Parms = np.zeros((NSteps, 24))  # 2D array of input parameters - for multiple voxels
for i in range(NSteps):
    Parms[i] = ParmLocal  # most of the parameters are the same in all voxels
    Parms[i, 4] = 50.0 + 30.0 * i / (NSteps - 1)  # the viewing angle varies from 50 to 80 degrees along the LOS

# parameters of the electron distribution function
Emin = 0.1  # E_min, MeV
Emax = 10.  # E_max, MeV
delta = 4.  # delta_1
n_b = 1e6  # n_b - nonthermal electron density, cm^{-3}
mu_c = np.cos(np.pi * 70. / 180.)  # loss-cone boundary
dmu_c = 0.2  # Delta_mu

E_arr = np.logspace(np.log10(Emin), np.log10(Emax), N_E)  # energy grid (log spaced)
mu_arr = np.linspace(-1., 1., N_mu)  # pitch-angle grid

f0 = np.zeros((N_mu, N_E))  # 2D distribution function array - for a single voxel

# computing the distribution function (equivalent to PLW & GLC)
A = n_b / (2. * np.pi) * (delta - 1.) / (Emin ** (1. - delta) - Emax ** (1. - delta))
B = 0.5 / (mu_c + dmu_c * np.sqrt(np.pi) / 2 * math.erf((1. - mu_c) / dmu_c))
for i in range(N_E):
    for j in range(N_mu):
        amu = abs(mu_arr[j])
        if amu < mu_c:
            mu_factor = 1.
        else:
            mu_factor = np.exp(-((amu - mu_c) / dmu_c) ** 2)
        f0[j, i] = A * B * E_arr[i] ** (-delta) * mu_factor

f_arr = np.zeros((NSteps, N_mu, N_E))  # 3D distribution function array - for multiple voxels
for k in range(NSteps):
    f_arr[k] = f0  # electron distribution function is the same in all voxels

# calculating the emission for analytical distribution (array -> on)
RL = pygs.get_mw(libname, Lparms, Rparms, Parms, E_arr, mu_arr, f_arr)

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
