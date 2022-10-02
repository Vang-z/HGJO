%% Add path
addpath(genpath('CEC/'));
addpath(genpath('Utils/'));

%% Initialize params
clc;

population_size = 60;
max_iteration = 1000;
N_function = 12;
N_algorithm = 1;
runtimes = 3;

for i_func = 1:N_function
    [lb, ub, dim, f_obj, f_name] = CEC22(i_func);

    for i_algorithm = 1:N_algorithm
        switch i_algorithm
            case 1
                % 2022
                algorithm = 'HGJO';
        end

        fprintf('================================== test function: %s ==================================\n', f_name);
        for i = 1:runtimes
            tic;
            [fitness, Convergence] = feval((algorithm), population_size, dim, lb, ub, max_iteration, f_obj);
            runtime = toc;
            T = clock;
            fprintf('[%s-%2s-%2s %2s:%2s:%2s]: %4s, runtimes: %2s, run time: %5fs, Fiteness: %f\n', ...
                    num2str(T(1)), num2str(T(2)), num2str(T(3)), num2str(T(4)), num2str(T(5)), num2str(floor(T(6))), ...
                    algorithm, num2str(i), runtime, fitness(end));
            Indicator(i_func).(algorithm)(i).fitness = fitness;
            Indicator(i_func).(algorithm)(i).Convergence = Convergence;
            clear fitness Convergence;
        end
        fitness = cell2mat({Indicator(i_func).(algorithm).fitness}');
        choosen_In = median(fitness, 2);
        Table((i_func - 1) * 5 + 1, i_algorithm) = max(choosen_In);
        Table((i_func - 1) * 5 + 2, i_algorithm) = min(choosen_In);
        Table((i_func - 1) * 5 + 3, i_algorithm) = mean(choosen_In);
        Table((i_func - 1) * 5 + 4, i_algorithm) = median(choosen_In);
        Table((i_func - 1) * 5 + 5, i_algorithm) = std(choosen_In);
        median_index = find(choosen_In == median(choosen_In));
        choose_fitness(i_func).(algorithm) = fitness(median_index(1), :);
        Convergence = cell2mat({Indicator(i_func).(algorithm).Convergence}');
        choose_Convergence(i_func).(algorithm) = Convergence(median_index(1), :);
        clear fitness Convergence;
    end


 %%  draw Convergence plot
     figure;
     semilogy(choose_Convergence(i_func).HGJO, 'linewidth', 1.5, 'color', '#d02725');
     grid on;

     title(f_name)
     xlabel('Iteration');
     ylabel('Best score obtained so far');
     axis tight;
     grid on;
     box on;
     legend('HGJO');

 %%  draw boxplot
     figure;
     HGJO_fitness = mean(cell2mat({Indicator(i_func).HGJO.fitness}'), 2);
     boxplot([HGJO_fitness], ...
             'Notch', 'on', 'Labels', {'HGJO'})
     title(f_name)
     xlabel('Algorithms')
     ylabel('Fitness Value')
     grid on;
end

save('CEC_Indicator', 'Indicator');

%exportgraphics(gca, 'F1_Convergence.png', 'Resolution', 96);
