function [intensity, Iout, prob, Convergence] = MT(Iin, level, N, Epochs, MA, f_func)

%% Initialize parameters
    % number of thresholds
    N_PAR = level  - 1;
    dim = N_PAR;
    % Calculate the total number of pixels in the image (ROW X COL)
    Nt = size(Iin, 1) * size(Iin, 2);
    % Lmax niveles de color a segmentar 1 - 256
    % 256 different maximum levels are considered in an image (i.e., 0 to 255)
    Lmax = 256;

    if size(Iin, 3) == 1
        % Imagen en escala de grises
        ub = ones(1, dim) * Lmax;
        lb = ones(1, dim);
    elseif size(Iin, 3) == 3
        % Imagen RGB
        ub = ones(1, dim) * Lmax;
        lb = ones(1, dim);
    end

    if size(Iin, 3) == 1
        % grayscale image
        [n_countR, x_valueR] = imhist(Iin(:, :, 1));
    elseif size(Iin, 3) == 3
        % RGB image
        % Calculate the color intensity histogram for each channel
        [n_countR, x_valueR] = imhist(Iin(:, :, 1));
        [n_countG, x_valueG] = imhist(Iin(:, :, 2));
        [n_countB, x_valueB] = imhist(Iin(:, :, 3));
    end

    % Calculate the probability distribution of color intensity
    for i = 1:Lmax
        if size(Iin, 3) == 1
            % grayscale image
            probR(i) = n_countR(i) / Nt;
        elseif size(Iin, 3) == 3
            % RGB image
            probR(i) = n_countR(i) / Nt;
            probG(i) = n_countG(i) / Nt;
            probB(i) = n_countB(i) / Nt;
        end
    end

    if size(Iin, 3) == 1
        % Initialize the population
        xR = chaotic(N, dim, ub, lb);
        for si=1:N
           xR(si, :) = sort(xR(si, :));
        end
        % calculate the fitness of population
        fitR = feval(f_func, N, level, xR, probR);
    elseif size(Iin,3) == 3
        % Initialize the population of the three channels
        xR = chaotic(N, dim, ub, lb);
        xG = chaotic(N, dim, ub, lb);
        xB = chaotic(N, dim, ub, lb);
        for si=1:N
           xR(si, :) = sort(xR(si, :));
           xG(si, :) = sort(xG(si, :));
           xB(si, :) = sort(xB(si, :));
        end
        fitR = feval(f_func, N, level, xR, probR);
        fitG = feval(f_func, N, level, xG, probG);
        fitB = feval(f_func, N, level, xB, probB);
    end

%% Iteration
    if size(Iin, 3) == 1
        [xR, fitR, Convergence] = feval(MA, level, N, dim, lb, ub, Epochs, xR, fitR, probR, f_func);
        [best_fitnessR, indexR] = max(fitR);
        best_xR = xR(indexR, :);
    elseif size(Iin, 3) == 3
        [xR, fitR, ConvergenceR] = feval(MA, level, N, dim, lb, ub, Epochs, xR, fitR, probR, f_func);
        [xG, fitG, ConvergenceG] = feval(MA, level, N, dim, lb, ub, Epochs, xG, fitG, probG, f_func);
        [xB, fitB, ConvergenceB] = feval(MA, level, N, dim, lb, ub, Epochs, xB, fitB, probB, f_func);

        [best_fitnessR, indexR] = max(fitR);
        [best_fitnessG, indexG] = max(fitG);
        [best_fitnessB, indexB] = max(fitB);
        best_xR = xR(indexR, :);
        best_xG = xR(indexG, :);
        best_xB = xR(indexB, :);
        Convergence = [ConvergenceR; ConvergenceG; ConvergenceB];
    end
%% Return the results
    if size(Iin, 3) == 1
        intensity = best_xR;
        Iout = img2gray(Iin, intensity);
        prob = [probR];
    elseif size(Iin, 3) == 3
        intensity = [best_xR; best_xG; best_xB];
        Iout = img2rgb(Iin, best_xR, best_xG, best_xB);
        prob = [probR; probG; probB];
    end
end
