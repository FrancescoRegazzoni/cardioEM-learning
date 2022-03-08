function [tt, V_log, p_log] = windkessel_model(parameters, rhs_LV, x0_LV)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Solve a windkessel model coupled with a left ventricle (LV) model. See [1] for more details.
%
% Inputs:
%    paramaters: Stucture containing the parameters of the simulation.
%    rhs_LV:     Right-hand side of the ODE system of LV dynamics. Inputs are the current
%                state, the current LV pressure (mmHg) and current time (s). Output is the
%                state rate of change.
%    x0_LV:      Initial state of the ODE system of LV dynamics. The first state coincides 
%                with the LV volume.
%
% Outputs:
%    tt:         Vector of time instants [s].
%    V_log:      Vector of LV volume [mL].
%    p_log:      Vector of PV pressures [mmHg].
%
% References:
%    [1] https://doi.org/10.1016/j.cma.2020.113268)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Windkessel model parameters
    AV_opening_pressure = parameters.AV_opening_pressure; % [Pa]
    MV_opening_pressure = parameters.MV_opening_pressure; % [Pa]
    tol                 = parameters.tol;                 % [mL/s]
    C                   = parameters.C;                   % [m^3 / Pa]
    R                   = parameters.R;                   % [Pa * s / m^3]

    %% Heartbeat parameters
    THB                 = parameters.THB;                 % [s]
    dt                  = parameters.dt;                  % [s]
    T_end               = parameters.T_end;               % [s]
    p_start             = parameters.p_start;             % [mmHg]
    
    %% Initialization
    mmHg_to_Pa = 133.322;
    AV_opening_pressure = AV_opening_pressure / mmHg_to_Pa; % [Pa] -> [mmHg]
    MV_opening_pressure = MV_opening_pressure / mmHg_to_Pa; % [Pa] -> [mmHg]
    C = C * mmHg_to_Pa * 1e6; % [m^3 / Pa] -> [mL / mmHg]
    R = R / mmHg_to_Pa / 1e6; % [Pa * s / m^3] -> [mmHg * s / mL]
    
    tt = 0:dt:T_end;
    phase = 1;
    
    x = x0_LV;
    V = x(1);
    p = p_start;
    
    p_log = zeros(length(tt), 1);
    V_log = zeros(length(tt), 1);
    p_log(1) = p;
    V_log(1) = V;
    
    t_last_change = 0;
    t_min_delay = 0.01; % [s]
    
    opts = optimoptions('fsolve','Algorithm','levenberg-marquardt', 'StepTolerance', 1e-6, 'FunctionTolerance', 1e-6, 'Display','off');
    
    %% Time loop
    h = waitbar(0,sprintf('Solving Windkessel model + LV ROM (%d heartbeats)...', round(T_end / THB)));
    for i = 2:length(tt)
        if phase == 1 || phase == 3 % isochoric phases
            y_new = fsolve(@(y) [ ...
                y(1:end-1) - x - dt * rhs_LV(y(1:end-1), y(end), tt(i)); % ANN model
                y(1) - V                                                 % volume constraint
                ], [x;p], opts);
            x = y_new(1:end-1);
            p = y_new(end);
            if phase == 1 && p > AV_opening_pressure && tt(i) >= t_last_change + t_min_delay
                phase = 2; % Aortic valve opens
                t_last_change = tt(i);
            elseif phase == 3 && p < MV_opening_pressure && tt(i) >= t_last_change + t_min_delay
                phase = 4; % Mitral valve opens
                t_last_change = tt(i);
            end
        elseif phase == 2 % ejection
            p = ((R * C) / (R * C + dt)) * p - (R / (R * C + dt)) * (V_log(i-1) - V_log(i-2)); % windkessel
            x = x + dt * rhs_LV(x, p, tt(i));
            if ((x(1) - V_log(i-1)) / dt > tol && tt(i) >= t_last_change + t_min_delay)
                phase = 3; % Aortic valve closes
                t_last_change = tt(i);
            end
        elseif phase == 4 % filling
            p = p + (p_start - p) / ((THB - mod(tt(i), THB))/dt); % linear ramp
            x = x + dt * rhs_LV(x, p, tt(i));
            if (mod(tt(i) + dt, THB) < 1e-1)
                phase = 1; % Mitral valve closes
                t_last_change = tt(i);
            end
        end
        V = x(1);
        p_log(i) = p;
        V_log(i) = V;
        waitbar(i / length(tt))
    end
    close(h)
end