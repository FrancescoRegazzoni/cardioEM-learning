clear

THB = 0.8; % [s]

use_cases{1}.problem = 'problems/EM_one_param.ini';
use_cases{1}.folder_name = 'simulation_one_param_';
use_cases{1}.samples_number = 35;

use_cases{2}.problem = 'problems/EM_all_params.ini';
use_cases{2}.folder_name = 'simulation_all_params_';
use_cases{2}.samples_number = 45;

for idx = 1:length(use_cases)
    problem = problem_get('app_cardioEM-learning', use_cases{idx}.problem);

    % Loop over the simulations contained in the 'data' folder.
    for i = 1:use_cases{idx}.samples_number
        fprintf('parsing simulation %d of %d...\n', i, use_cases{idx}.samples_number)
        current_folder = ['data/' use_cases{idx}.folder_name sprintf('%02d', i)];
        T = readtable([current_folder '/output.csv']);
        params = jsondecode(fileread([current_folder '/parameters.json']));

        tt = T.time';           % [s]
        uu = [cos(2*pi*tt/THB); % [-]
              sin(2*pi*tt/THB); % [-]
              T.pressure'];     % [mmHg]
 
        % The parameters of the EM models are seen as time-dependent inputs
        % that are however constant in time.
        switch use_cases{idx}.problem 
            case 'problems/EM_one_param.ini'
                uu = [uu;
                      ones(1,length(tt)) * params.a_XB];
            case 'problems/EM_all_params.ini'
                uu = [uu;
                      ones(1,length(tt)) * params.a_XB;
                      ones(1,length(tt)) * params.sigma_f;
                      ones(1,length(tt)) * params.alpha_fibers;
                      ones(1,length(tt)) * params.C_mechanics];
            otherwise
                error('non recognized problem type')
        end
        yy =  T.volume'; % [mV]

        dataset{i}.tt = tt;
        dataset{i}.uu = uu;
        dataset{i}.yy = yy;
    end

    dataset_save(problem, dataset, 'samples')
    fprintf('Generated dataset for problem %s\n', use_cases{idx}.problem)
end