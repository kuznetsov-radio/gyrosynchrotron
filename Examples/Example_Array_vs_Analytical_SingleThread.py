import numpy as np
import matplotlib.pyplot as plt
import GScodes
import math
import platform,sys,os
import warnings
homedir = os.path.expanduser("~")
if platform.system() == 'Darwin':
    if sys.version_info[0] < 3:
        sys.path.append(homedir + '/Dropbox/bc_python/gs/py27gs_mac/')
    if sys.version_info[0] >= 3:
        if sys.version_info[1] == 5:
            sys.path.append(homedir + '/Dropbox/bc_python/gs/py35gs_mac/')
        if sys.version_info[1] == 6:
            sys.path.append(homedir + '/Dropbox/bc_python/gs/py36gs_mac/')
            if sys.version_info[2] != 1:
                warnings.warn("GS codes only run in Python 3.6.1")
if platform.system() == 'Linux':
    sys.path.append(homedir + '/Dropbox/bc_python/gs/py36gs_linux/')

import MWTransfer

libname = '../Binaries/MWTransferArr.so'
GET_MW=GScodes.initGET_MW(libname)  # load the library
isotropic = True

# ============== Analytical Version ====================
Nf = 100  # number of frequencies
NSteps = 10  # number of nodes along the line-of-sight
N_E = 15  # number of energy nodes
N_mu = 15  # number of pitch-angle nodes

Lparms = np.zeros(11, dtype='int32')  # array of dimensions etc.
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
ParmLocal[4] = 80.  # theta - viewing angle, deg
ParmLocal[6] = 3  # distribution over energy (PLW is chosen)
ParmLocal[7] = 1e6  # n_b - nonthermal electron density, cm^{-3}
ParmLocal[9] = 0.1  # E_min, MeV
ParmLocal[10] = 10.0  # E_max, MeV
ParmLocal[12] = 4.0  # \delta_1
if isotropic:
    ParmLocal[14] = 0  # distribution over pitch-angle (isotropic is chosen)
else:
    ParmLocal[14] = 3  # distribution over pitch-angle (GLC is chosen)
ParmLocal[15] = 70  # loss-cone boundary, degrees
ParmLocal[16] = 0.2  # \Delta\mu


Parms = np.zeros((24, NSteps), dtype='double', order='F')  # 2D array of input parameters - for multiple voxels
for i in range(NSteps):
    Parms[:, i] = ParmLocal  # most of the parameters are the same in all voxels
    #Parms[i, 4] = 50.0 + 30.0 * i / (NSteps - 1)  # the viewing angle varies from 50 to 80 degrees along the LOS

RL = np.zeros((7, Nf), dtype='double', order='F')  # input/output array
dummy = np.array(0, dtype='double')

# calculating the emission for analytical distribution (array -> off),
# the unused parameters can be set to any value
res = GET_MW(Lparms, Rparms, Parms, dummy, dummy, dummy, RL)

# plot the results
f_anal = RL[0]  # emission frequency, GHz
I_L = RL[5]  # left-polarized emission intensity (as observed from the Earth), sfu
I_R = RL[6]  # right-polarized emission intensity (as observed from the Earth), sfu
flux_anal = I_L + I_R
pol_anal = (I_L - I_R) / (I_L + I_R)

# ========================= Array Version ==============================
Lparms = np.zeros(11, dtype='int32')  # array of dimensions etc.
Lparms[0] = NSteps
Lparms[1] = Nf
Lparms[2] = N_E
Lparms[3] = N_mu
if isotropic:
    Lparms[9] = 1 # force isotropic

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
ParmLocal[4] = 80.  # theta - viewing angle, deg

Parms = np.zeros((24, NSteps), dtype='double', order='F')  # 2D array of input parameters - for multiple voxels
for i in range(NSteps):
    Parms[:, i] = ParmLocal  # most of the parameters are the same in all voxels
    #Parms[i, 4] = 50.0 + 30.0 * i / (NSteps - 1)  # the viewing angle varies from 50 to 80 degrees along the LOS

# parameters of the electron distribution function
Emin = 0.1  # E_min, MeV
Emax = 10.  # E_max, MeV
delta = 4.  # delta_1
n_b = 1e6  # n_b - nonthermal electron density, cm^{-3}
mu_c = np.cos(np.pi * 70. / 180.)  # loss-cone boundary
dmu_c = 0.2  # Delta_mu

E_arr = np.logspace(np.log10(Emin), np.log10(Emax), N_E)  # energy grid (log spaced)
mu_arr = np.linspace(-1., 1., N_mu)  # pitch-angle grid

f0 = np.zeros((N_E, N_mu))  # 2D distribution function array - for a single voxel

# computing the distribution function (equivalent to PLW & GLC)
A = n_b/(2.0*np.pi)*(delta-1.0)/(Emin**(1.0-delta)-Emax**(1.0-delta))
B = 0.5/(mu_c+dmu_c*np.sqrt(np.pi)/2*math.erf((1.0-mu_c)/dmu_c))
for i in range(N_E):
    for j in range(N_mu):
        amu = abs(mu_arr[j])
        f0[i, j] = A*B*E_arr[i]**(-delta)*(1.0 if amu < mu_c else np.exp(-((amu-mu_c)/dmu_c)**2))

f_arr = np.zeros((N_E, N_mu, NSteps), dtype='double', order='F')  # 3D distribution function array - for multiple voxels
for k in range(NSteps):
    f_arr[:, :, k] = f0  # electron distribution function is the same in all voxels

RL = np.zeros((7, Nf), dtype='double', order='F')  # input/output array

# calculating the emission for analytical distribution (array -> on)
res = GET_MW(Lparms, Rparms, Parms, E_arr, mu_arr, f_arr, RL)

# plot the results
f_arr = RL[0]  # emission frequency, GHz
I_L = RL[5]  # left-polarized emission intensity (as observed from the Earth), sfu
I_R = RL[6]  # right-polarized emission intensity (as observed from the Earth), sfu
flux_arr = I_L + I_R
pol_arr = (I_L - I_R) / (I_L + I_R)

# ================= Former Analytical Version =================
ParmLocal = 29 * [0]  # the array of the volume element parameters
ParmLocal[0] = 1e20  # Area, cm^2
ParmLocal[1] = 1e10  # Depth, cm (will be changed later!)
ParmLocal[2] = 3e7  # T_0, K
ParmLocal[3] = 0.05  # \eps (not used in this example)
ParmLocal[4] = 4.0  # \kappa (not used in this example)
ParmLocal[5] = 16  # number of integration nodes
ParmLocal[6] = 0.1  # E_min, MeV
ParmLocal[7] = 10.0  # E_max, MeV
ParmLocal[8] = 1.0  # E_break, MeV (not used in this example)
ParmLocal[9] = 4.0  # \delta_1
ParmLocal[10] = 6.0  # \delta_2 (not used in this example)
ParmLocal[11] = 3e9  # n_0 - thermal electron density, cm^{-3}
ParmLocal[12] = 1e6  # n_b - nonthermal electron density, cm^{-3}
ParmLocal[13] = 180  # B - magnetic field, G
ParmLocal[14] = 80.  # theta - the viewing angle, degrees
ParmLocal[15] = 1e9  # starting frequency to calculate spectrum, Hz
ParmLocal[16] = 0.02  # logarithmic step in frequency
ParmLocal[17] = 3  # distribution over energy (PLW is chosen)
ParmLocal[18] = Nf  # number of frequencies (specified above)
if isotropic:
    ParmLocal[19] = 1  # distribution over pitch-angle (isotropic is chosen)
else:
    ParmLocal[19] = 3  # distribution over pitch-angle (GLC is chosen)
ParmLocal[20] = 70  # loss-cone boundary, degrees
ParmLocal[21] = 0  # beam direction (degrees) in GAU and SGA (not used in this example)
ParmLocal[22] = 0.2  # \Delta\mu
ParmLocal[23] = 1  # a_4 in SGA (not used in this example)
ParmLocal[25] = 12  # f^C_cr
ParmLocal[26] = 12  # f^WH_cr
ParmLocal[27] = 1  # matching on
ParmLocal[28] = 1  # Q-optimization on

if sys.version_info[0] < 3:
    for i, p in enumerate(ParmLocal):
        if type(p) == int:
            ParmLocal[i] = long(p)

Parms0 = NSteps * [None]  # the array of input parameters
for i in range(NSteps):
    Parms0[i] = ParmLocal[:]  # most of the parameters are the same in all voxels
    Parms0[i][1] = ParmLocal[1] / NSteps  # But the length of an elementary interval should equal the source depth divided by nsteps
    #Parms0[i][14] = 50.0 + 30.0 * i / (NSteps - 1)  # the viewing angle varies from 50 to 80 degrees along the LOS

res = MWTransfer.GET_MW(Parms0)
f_anal_former = np.array(res[0])  # emission frequency (GHz)
I_L = np.array(res[1])  # left-hand polarized emission intensity, sfu (as observed at the Earth)
I_R = np.array(res[2])  # right-hand polarized emission intensity, sfu (as observed at the Earth)

flux_anal_former = I_L + I_R
pol_anal_former = (I_L - I_R) / (I_L + I_R)

# Plot the results
plt.figure(figsize=(12, 10))
# first plot, total intensity
ax = plt.subplot(221)
ax.plot(f_anal_former, flux_anal_former, color='k', label='Former PyAnalytical')
ax.plot(f_anal, flux_anal, color='b', label='New PyAnalytical')
ax.plot(f_arr, flux_arr, color='r', label='New PyArray')
ax.set_xlabel('Frequency (GHz)')
ax.set_ylabel('Flux (sfu)')
ax.set_title('Total Intensity')
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend()

# second plot, circular polarization degree
ax = plt.subplot(222)
ax.plot(f_anal_former, pol_anal_former, color='k', label='Former PyAnalytical')
ax.plot(f_anal, pol_anal, color='b', label='New PyAnalytical')
ax.plot(f_arr, pol_arr, color='r', label='New PyArray')
ax.set_xlabel('Frequency (GHz)')
ax.set_ylabel('Polarization Degree')
ax.set_title('Circular Polarization Degree')
ax.set_xscale('log')
ax.set_yscale('linear')
ax.set_ylim([-1, 1])
ax.legend()

# 3rd plot, total intensity ratios
ax = plt.subplot(223)
ax.plot(f_anal, np.ones_like(f_anal), color='k', label='Former PyAnalytical')
ax.plot(f_anal_former, flux_anal / flux_anal_former, color='b', label='New PyAnalytical')
ax.plot(f_arr, flux_arr / flux_anal_former, color='r', label='New PyArray')
ax.set_xlabel('Frequency (GHz)')
ax.set_ylabel('Ratio')
ax.set_title('Total Intensity Ratio')
ax.set_xscale('log')
ax.set_yscale('linear')
ax.set_ylim([0., 2.])
ax.legend()

# 4th plot, polarization degree ratios
ax = plt.subplot(224)
ax.plot(f_anal, np.ones_like(f_anal), color='k', label='Former PyAnalytical')
ax.plot(f_anal_former, pol_anal / pol_anal_former, color='b', label='New PyAnalytical')
ax.plot(f_arr, pol_arr / pol_anal_former, color='r', label='New PyArray')
ax.set_xlabel('Frequency (GHz)')
ax.set_ylabel('Ratio')
ax.set_title('Polarization Degree Ratio')
ax.set_xscale('log')
ax.set_yscale('linear')
ax.set_ylim([0., 2.])
ax.legend()

plt.tight_layout()
plt.show()
