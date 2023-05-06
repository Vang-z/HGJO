function [fitness, Convergence] = HGJO(N, dim, lb, ub, Epochs, f_obj)

%% Initialize parameters
    x = lb + rand(N, dim) .* (ub - lb);
    fitness = f_obj(x);
    % the OBL strategy
    x_obl = lb + ub - x;
    fitness_obl = f_obj(x_obl);
    fitness = [fitness, fitness_obl];
    x = [x; x_obl];
    [fitness, index] = sort(fitness);
    fitness = fitness(1:N);
    x = x(index(1:N), :);

    Male_jackal = zeros(1, dim);
    Male_jackal_fit = inf;

    Female_jackal = zeros(1, dim);
    Female_jackal_fit = inf;

    Convergence = zeros(1, Epochs);

%% Iteration
    for epoch = 1:Epochs
        rl = 0.05 * levy(N, dim, 1.5);
        % cauchy distribution and the Escape energy
        rc = cauchy([N, dim], 0, 0.5);
        E1 = 2 * (1 - epoch / Epochs)^(pi * epoch / Epochs);

        for i = 1:N
            if fitness(i) < Male_jackal_fit
                Male_jackal_fit = fitness(i);
                Male_jackal = x(i, :);
            elseif fitness(i) < Female_jackal_fit
                Female_jackal_fit = fitness(i);
                Female_jackal = x(i, :);
            end
        end

        for i = 1:N
            E0 = 2 * rand(1, dim) - 1;
            E = E1 * E0;
            if mean(abs(E)) > 1
                if rand > 0.5
                    Y1 = x(i, :) - E .* (Male_jackal - rc(i, :) .* x(i, :));
                    Y2 = x(i, :) - E .* (Female_jackal - rc(i, :) .* x(i, :));
                else
                    Y1 = x(i, :) - E .* (Male_jackal - rc(i, :) .* mean(x));
                    Y2 = x(i, :) - E .* (Female_jackal - rc(i, :) .* mean(x));
                end
            elseif mean(abs(E)) > 0.5
                if rand > 0.5
                    Y1 = x(i, :) - E .* (Male_jackal - rc(i, :) .* x(i, :));
                    Y2 = x(i, :) - E .* (Female_jackal - rc(i, :) .* x(i, :));
                else
                    Y1 = Male_jackal - E .* (rl(i, :) .* Male_jackal - x(i, :));
                    Y2 = Female_jackal - E .* (rl(i, :) .* Female_jackal - x(i, :));
                end
            else
                if rand > 0.5
                    Y1 = Male_jackal - E .* (rl(i, :) .* Male_jackal - x(i, :));
                    Y2 = Female_jackal - E .* (rl(i, :) .* Female_jackal - x(i, :));
                else
                    Y1 = Male_jackal - E .* (rl(i, :) .* Male_jackal - mean(x));
                    Y2 = Female_jackal - E .* (rl(i, :) .* Female_jackal - mean(x));
                end
            end
            Y = (Y1 + Y2) / 2;
            Flag4lb = Y < lb;
            Flag4ub = Y > ub;
            Y = Y .* (~(Flag4ub + Flag4lb)) + ub .* Flag4ub + lb .* Flag4lb;

            fit_y = f_obj(Y);
            if fit_y < fitness(i)
                x(i, :) = Y;
                fitness(i) = fit_y;
            end

            % combine helpers
            helper_idx = randperm(N, 3);
            helper1 = x(helper_idx(1), :);
            helper2 = x(helper_idx(2), :);
            helper3 = x(helper_idx(3), :);
            offspring = helper1 + E .* (helper2 - helper3);
            Flag4lb = offspring < lb;
            Flag4ub = offspring > ub;
            offspring = offspring .* (~(Flag4ub + Flag4lb)) + ub .* Flag4ub + lb .* Flag4lb;
            fit = f_obj(offspring);
            if fit < fitness(i)
                x(i, :) = offspring;
                fitness(i) = fit;
            end
        end

        % helpers
        helpers = x - rand * (x(randperm(N), :) - x(randperm(N), :));
        Flag4lb = helpers < lb;
        Flag4ub = helpers > ub;
        helpers = helpers .* (~(Flag4ub + Flag4lb)) + ub .* Flag4ub + lb .* Flag4lb;
        helpers_fit = f_obj(helpers);
        flag4helpers = (helpers_fit < fitness);
        u = repmat(flag4helpers', 1, dim);
        x = u .* helpers + ~u .* x;
        fitness = flag4helpers .* helpers_fit + ~flag4helpers .* fitness;

        Convergence(epoch) = min(fitness);
    end
end
