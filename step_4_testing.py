#!/usr/bin/env python 
use_case = 'one_param' # 'one_param' or 'all_params'

#%%
from pyModelLearning import ANNmodel
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
from windkessel_model import windkessel_model

# %% Windkessel model parameters
parameters = dict()
parameters['AV_opening_pressure'] = 11000  # [Pa]
parameters['MV_opening_pressure'] = 667    # [Pa]
parameters['tol']                 = 1e-5   # [mL/s]
parameters['C']                   = 4.5e-9 # [m^3 / Pa]
parameters['R']                   = 5.0e7  # [Pa * s / m^3]

# %% Heartbeat parameters
parameters['THB']                 = 0.8    # [s]
parameters['T_end']               = 4      # [s]
parameters['p_start']             = 16.8   # [mmHg]
parameters['dt']                  = 1e-3   # [s]

# %% LV parameters
aXB                               = 160    # [MPa]
sigma_f                           = 76.43  # [mm^2 / s]
alpha_fibers                      = 60     # [deg]
C_mechanics                       = 880    # [Pa]

# %% Loading ANN model
if use_case == 'one_param':
    model_dir = 'app_cardioEM-learning/networks_EM_one_param/ROM_int_N2_hlayF8_dof74_ntrain24/compact.mat'
    param_vec = np.array([aXB])
elif use_case == 'all_params':
    model_dir = 'app_cardioEM-learning/networks_EM_all_params/ROM_int_N1_hlayF12_dof121_ntrain32/compact.mat'
    param_vec = np.array([aXB, sigma_f, alpha_fibers, C_mechanics])

ANNmod = ANNmodel(model_dir)
THB = parameters['THB']
rhs_LV = lambda x, p, t: ANNmod.rhs(x, np.concatenate([np.array([np.cos(2*np.pi*t/THB), np.sin(2*np.pi*t/THB), p]), param_vec]))
x0_LV = ANNmod.initial_state

print('solving Windkessel model + LV ROM...', end = '')
tt, V_log, p_log = windkessel_model(parameters, rhs_LV, x0_LV)
print(' done!')

#%% Reference solution (full-order model)
T = pd.read_csv('data/simulation_windkessel/output.csv')

# %% compute errors
mmHg_to_Pa = 133.322

# extract last cycle
idxs_ROM = tt > 3.2
idxs_FOM = T.time > 3.2
t_ROM = tt[idxs_ROM]
t_FOM = T.time[idxs_FOM]
p_ROM = p_log[idxs_ROM]
p_FOM = T.pressure[idxs_FOM]
V_ROM = V_log[idxs_ROM]
V_FOM = T.volume[idxs_FOM]

p_FOM_interp = interpolate.interp1d(t_FOM, p_FOM)(t_ROM)
V_FOM_interp = interpolate.interp1d(t_FOM, V_FOM)(t_ROM)

p_max_ROM = max(p_ROM)
p_max_FOM = max(p_FOM)
p_min_ROM = min(p_ROM)
p_min_FOM = min(p_FOM)
V_max_ROM = max(V_ROM)
V_max_FOM = max(V_FOM)
V_min_ROM = min(V_ROM)
V_min_FOM = min(V_FOM)

print('')
print('Computed biomarkers:')
print('p_max (FOM): %f mmHg' % p_max_FOM)
print('p_max (ROM): %f mmHg' % p_max_ROM)
print('p_min (FOM): %f mmHg' % p_min_FOM)
print('p_min (ROM): %f mmHg' % p_min_ROM)
print('V_max (FOM): %f mV'   % V_max_FOM)
print('V_max (ROM): %f mV'   % V_max_ROM)
print('V_min (FOM): %f mV'   % V_min_FOM)
print('V_min (ROM): %f mV'   % V_min_ROM)

print('')
print('Errors on transients (relative L^2 error):')
print('pressure: %1.2e' % (np.linalg.norm(p_ROM-p_FOM_interp) / np.linalg.norm(p_FOM_interp)))
print('volume:   %1.2e' % (np.linalg.norm(V_ROM-V_FOM_interp) / np.linalg.norm(V_FOM_interp)))

print('')
print('Errors on biomarkers:')
print('p_min: relative %1.2e  absolute %f mmHg' % (abs(p_min_ROM-p_min_FOM)/p_min_FOM, abs(p_min_ROM-p_min_FOM)))
print('p_max: relative %1.2e  absolute %f mmHg' % (abs(p_max_ROM-p_max_FOM)/p_max_FOM, abs(p_max_ROM-p_max_FOM)))
print('V_min: relative %1.2e  absolute %f mL' % (abs(V_min_ROM-V_min_FOM)/V_min_FOM, abs(V_min_ROM-V_min_FOM)))
print('V_max: relative %1.2e  absolute %f mL' % (abs(V_max_ROM-V_max_FOM)/V_max_FOM, abs(V_max_ROM-V_max_FOM)))

# %% postprocessing

fig = plt.figure(figsize=(6,3))

ax_V = plt.subplot2grid((2, 2), (0, 0))
ax_p = plt.subplot2grid((2, 2), (1, 0))
ax_pV = plt.subplot2grid((2, 2), (0, 1), rowspan=2)

ax_V.plot(T.time, T.volume)
ax_V.plot(tt, V_log, '--')
ax_V.set_xlabel('t [s]')
ax_V.set_ylabel('V [mL]')
ax_V.set_xlim([0, parameters['T_end']])
ax_V.set_ylim([60, 150])

ax_p.plot(T.time, T.pressure)
ax_p.plot(tt, p_log, '--')
ax_p.set_xlabel('t [s]')
ax_p.set_ylabel('p [mmHg]')
ax_p.set_xlim([0, parameters['T_end']])
ax_p.set_ylim([0, 120])

ax_pV.plot(T.volume, T.pressure, label = 'FOM')
ax_pV.plot(V_log, p_log, '--', label = 'ROM')
ax_pV.set_xlabel('V [mL]')
ax_pV.set_ylabel('p [mmHg]')
ax_pV.legend()
ax_pV.set_xlim([60, 150])
ax_pV.set_ylim([0, 120])

fig.tight_layout()

plt.show()