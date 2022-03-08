import numpy as np
from scipy.optimize import root
from scipy import interpolate

def windkessel_model(parameters, rhs_LV, x0_LV):
###############################################################################################
#
# Solve a windkessel model coupled with a left ventricle (LV) model. See [1] for more details.
#
# Inputs:
#    paramaters: Dictionary containing the parameters of the simulation.
#    rhs_LV:     Right-hand side of the ODE system of LV dynamics. Inputs are the current
#                state, the current LV pressure (mmHg) and current time (s). Output is the
#                state rate of change.
#    x0_LV:      Initial state of the ODE system of LV dynamics. The first state coincides 
#                with the LV volume.
#
# Outputs:
#    tt:         Vector of time instants [s].
#    V_log:      Vector of LV volume [mL].
#    p_log:      Vector of PV pressures [mmHg].
#
# References:
#    [1] https://doi.org/10.1016/j.cma.2020.113268)
#
###############################################################################################

    # %% Windkessel model parameters
    AV_opening_pressure = parameters['AV_opening_pressure'] # [Pa]
    MV_opening_pressure = parameters['MV_opening_pressure'] # [Pa]
    tol                 = parameters['tol']                 # [mL/s]
    C                   = parameters['C']                   # [m^3 / Pa]
    R                   = parameters['R']                   # [Pa * s / m^3]

    # %% Heartbeat parameters
    THB                 = parameters['THB']                 # [s]
    dt                  = parameters['dt']                  # [s]
    T_end               = parameters['T_end']               # [s]
    p_start             = parameters['p_start']             # [mmHg]
    
    # %% Initialization
    mmHg_to_Pa = 133.322
    AV_opening_pressure = AV_opening_pressure / mmHg_to_Pa
    MV_opening_pressure = MV_opening_pressure / mmHg_to_Pa
    C = C * mmHg_to_Pa * 1e6
    R = R / mmHg_to_Pa / 1e6

    tt = np.arange(0, T_end + 1e-10*dt, dt)
    phase = 1

    x = x0_LV
    V = x[0]
    p = p_start

    p_log = np.zeros((len(tt)))
    V_log = np.zeros((len(tt)))
    p_ED = p
    p_log[0] = p
    V_log[0] = V
    
    t_last_change = 0
    t_min_delay = 0.01 # [s]

    # %% Time loop
    for i in range(1, len(tt)):
        if phase == 1 or phase == 3: # isochoric phases
            residual = lambda y: np.concatenate([
                    y[:-1] - x - dt * rhs_LV(y[:-1], y[-1], tt[i]), # ANN model
                    np.array([y[0] - V])                            # volume constraint
                ])
            res = root(residual, np.concatenate([x, np.array([p])]), method = 'lm')
            x = res.x[:-1]
            p = res.x[-1]
            if phase == 1 and p > AV_opening_pressure and tt[i] >= t_last_change + t_min_delay:
                phase = 2 # Aortic valve opens
                t_last_change = tt[i]
            elif phase == 3 and p < MV_opening_pressure and tt[i] >= t_last_change + t_min_delay:
                phase = 4 # Mitral valve opens
                t_last_change = tt[i]
        elif phase == 2: # ejection
            p = ((R * C) / (R * C + dt)) * p - (R / (R * C + dt)) * (V_log[i-1] - V_log[i-2]) # windkessel
            x = x + dt * rhs_LV(x, p, tt[i])
            if ((x[0] - V_log[i-1]) / dt > tol and tt[i] % THB > 0.01) and tt[i] >= t_last_change + t_min_delay:
                phase = 3 # Aortic valve closes
                t_last_change = tt[i]
        elif phase == 4: # filling
            p = p + (p_ED - p) / ((THB - tt[i] % THB)/dt) # linear ramp
            x = x + dt * rhs_LV(x, p, tt[i])
            if ((tt[i] + dt) % THB < 1e-1):
                phase = 1 # Mitral valve closes
                t_last_change = tt[i]
        V = x[0]
        p_log[i] = p
        V_log[i] = V

    return tt, V_log, p_log