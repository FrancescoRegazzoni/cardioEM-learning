clear
use_case = 'one_param'; % 'one_param' or 'all_params';

%% Windkessel model parameters
parameters.AV_opening_pressure = 11000;  % [Pa]
parameters.MV_opening_pressure = 667;    % [Pa]
parameters.tol                 = 1e-5;   % [mL/s]
parameters.C                   = 4.5e-9; % [m^3 / Pa]
parameters.R                   = 5.0e7;  % [Pa * s / m^3]

%% Heartbeat parameters
parameters.THB                 = 0.8;    % [s]
parameters.T_end               = 4;      % [s]
parameters.p_start             = 16.8;   % [mmHg]
parameters.dt                  = 1e-3;   % [s]

%% LV parameters
aXB                            = 160;    % [MPa]
sigma_f                        = 76.43;  % [mm^2 / s]
alpha_fibers                   = 60;     % [deg]
C_mechanics                    = 880;    % [Pa]

%% Loading ANN model
switch use_case
    case 'one_param'
        problem_name = 'problems/EM_one_param.ini';
        model_dir = 'ROM_int_N2_hlayF8_dof74_ntrain24';
        param_vec = aXB;
    case 'all_params'
        problem_name = 'problems/EM_all_params.ini';
        model_dir = 'ROM_int_N1_hlayF12_dof121_ntrain32';
        param_vec = [aXB; sigma_f; alpha_fibers; C_mechanics];
end

problem = problem_get('app_cardioEM-learning', problem_name);
ANNmod = read_model_fromfile(problem, model_dir);
rhs_LV = @(x, p, t) ANNmod.f(x, [cos(2*pi*t/parameters.THB); sin(2*pi*t/parameters.THB); p; param_vec]);
x0_LV = ANNmod.x0;
fprintf('solving Windkessel model + LV ROM...')
[tt, V_log, p_log] = windkessel_model(parameters, rhs_LV, x0_LV);
fprintf(' done!\n')

%% Reference solution (full-order model)
T = readtable('data/simulation_windkessel/output.csv');

%% compute errors
mmHg_to_Pa = 133.322;

% extract last cycle
idxs_ROM = tt > 3.2;
idxs_FOM = T.time > 3.2;
t_ROM = tt(idxs_ROM);
t_FOM = T.time(idxs_FOM);
p_ROM = p_log(idxs_ROM);
p_FOM = T.pressure(idxs_FOM);
V_ROM = V_log(idxs_ROM);
V_FOM = T.volume(idxs_FOM);

p_FOM_interp = interp1(t_FOM, p_FOM, t_ROM)';
V_FOM_interp = interp1(t_FOM, V_FOM, t_ROM)';

p_max_ROM = max(p_ROM);
p_max_FOM = max(p_FOM);
p_min_ROM = min(p_ROM);
p_min_FOM = min(p_FOM);
V_max_ROM = max(V_ROM);
V_max_FOM = max(V_FOM);
V_min_ROM = min(V_ROM);
V_min_FOM = min(V_FOM);

fprintf('\nComputed biomarkers:\n')
fprintf('p_max (FOM): %f mmHg\n' , p_max_FOM)
fprintf('p_max (ROM): %f mmHg\n' , p_max_ROM)
fprintf('p_min (FOM): %f mmHg\n' , p_min_FOM)
fprintf('p_min (ROM): %f mmHg\n' , p_min_ROM)
fprintf('V_max (FOM): %f mV\n'   , V_max_FOM)
fprintf('V_max (ROM): %f mV\n'   , V_max_ROM)
fprintf('V_min (FOM): %f mV\n'   , V_min_FOM)
fprintf('V_min (ROM): %f mV\n'   , V_min_ROM)

fprintf('\nErrors on transients (relative L^2 error):\n')
fprintf('pressure: %1.2e\n', norm(p_ROM-p_FOM_interp) / norm(p_FOM_interp))
fprintf('volume:   %1.2e\n', norm(V_ROM-V_FOM_interp) / norm(V_FOM_interp))

fprintf('\nErrors on biomarkers:\n')
fprintf('p_min:  relative %1.2e  absolute %f mmHg\n', abs(p_min_ROM-p_min_FOM)/p_min_FOM, abs(p_min_ROM-p_min_FOM))
fprintf('p_max:  relative %1.2e  absolute %f mmHg\n', abs(p_max_ROM-p_max_FOM)/p_max_FOM, abs(p_max_ROM-p_max_FOM))
fprintf('V_min:  relative %1.2e  absolute %f mL\n'  , abs(V_min_ROM-V_min_FOM)/V_min_FOM, abs(V_min_ROM-V_min_FOM))
fprintf('V_max:  relative %1.2e  absolute %f mL\n'  , abs(V_max_ROM-V_max_FOM)/V_max_FOM, abs(V_max_ROM-V_max_FOM))

%% post-processing
figure('Position',  [100, 100, 600, 300])

subplot(2,2,1)
plot(T.time, T.volume, 'linewidth', 2)
hold on
plot(tt, V_log, '--', 'linewidth', 2);
xlabel('t [s]')
ylabel('V [mL]')
xlim([0, parameters.T_end])
ylim([60, 150])

subplot(2,2,3)
plot(T.time, T.pressure, 'linewidth', 2)
hold on
plot(tt, p_log, '--', 'linewidth', 2);
xlabel('t [s]')
ylabel('p [mmHg]')
xlim([0, parameters.T_end])
ylim([0, 120])

subplot(2,2,[2 4])
plot(T.volume, T.pressure, 'linewidth', 2)
hold on
plot(V_log, p_log, '--', 'linewidth', 2);
xlabel('V [mL]')
ylabel('p [mmHg]')
xlim([60, 150])
ylim([0, 120])

legend('FOM', 'ROM')