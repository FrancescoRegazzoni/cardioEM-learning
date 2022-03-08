% select here the use case
use_case = 'one_param'; % 'one_param' or 'all_params';

%%
switch use_case
    case 'one_param'
        problem_name = 'problems/EM_one_param.ini';
    case 'all_params'
        problem_name = 'problems/EM_all_params.ini';
end

%%
problem = problem_get('app_cardioEM-learning', problem_name);

dataset_def.problem = problem;
dataset_def.type = 'file';
dataset_def.source = 'samples.mat;:';

dataset = dataset_get(dataset_def);

dataset_plot(dataset, problem, struct('normalized_u', 0, 'normalized_y', 0))
%%
figure()
for i = 1:length(dataset)
    plot(dataset{i}.yy(1,:), dataset{i}.uu(3,:), '-', 'linewidth', 2)
    hold on
end
xlabel('LV volume [mL]')
ylabel('LV pressure [mmHg]')