%% Add path
addpath(genpath('Multi-Thresholding/'));
addpath(genpath('Utils/'));

%% Initialize params
clear;
clc;
% Choose the fitness functions: otsu | kapur
f_func = 'otsu';
levels = [8, 16, 24, 32];
images = {'c0080.png', 'c0088.png', 'c0132.png', 'c0135.png', 'c0180.png', 'c0536.png', 'c1088.png', 'l0032.png', 'l0064.png', 'l0135.png', 'l0158.png', 'l0226.png', 'l0699.png', 'l0879.png', 'l1074.png', 'x0017.png'};
population_size = 60;
max_iteration = 200;
N_algorithm = 1;
runtimes = 31;

for i_img = 1:2
img = cell2mat(images(i_img));
for i_level = 1:length(levels)
level = levels(i_level);
Iin = imread(img);
set(0, 'defaultfigurecolor', 'w');

if size(Iin, 3) == 1
    error('This code only for rgb images, if you want test the gray image, pls modify it by your self!')
elseif size(Iin, 3) == 3
    figure();
    mesh3('channel', Iin(:, :, 1), 'theme', 'red');
    title('The red channel of the raw image');
    figure();
    mesh3('channel', Iin(:, :, 2), 'theme', 'green');
    title('The green channel of the raw image');
    figure();
    mesh3('channel', Iin(:, :, 3), 'theme', 'blue');
    title('The blue channel of the raw image');
end

%hist2d(Iin)
for i_algorithm = 1:N_algorithm
switch i_algorithm
    case 1
        algorithm = 'HGJO';
end

for i = 1:runtimes
    T = clock;
    fprintf('========================== [%s-%s-%s %s:%s:%s]: %s Running, runtimes: %d ==========================\n', ...
        num2str(T(1)), num2str(T(2)), num2str(T(3)), num2str(T(4)), num2str(T(5)), num2str(floor(T(6))), algorithm, i);
    tic;
    [intensity, Iout, prob, Convergence] = MT(Iin, level, population_size, max_iteration, algorithm, f_func);
    runtime = toc;
    fprintf('Run time: %12fs\n', runtime);

    psnr4mt = psnr(Iin, Iout);
    ssim4mt = ssim(Iin, Iout);
    [fsim, fsim_c] = FSIM(Iin, Iout);
    Indicator.(algorithm)(i).intensity = intensity;
    Indicator.(algorithm)(i).Iout = Iout;
    Indicator.(algorithm)(i).Convergence = Convergence;
    Indicator.(algorithm)(i).Fitness = Convergence(:, end);
    Indicator.(algorithm)(i).PSNR = psnr4mt;
    Indicator.(algorithm)(i).SSIM = ssim4mt;
    Indicator.(algorithm)(i).FSIM = fsim;
    Indicator.(algorithm)(i).FSIMc = fsim_c;

    disp('The best fitness: ');
    disp(['    Red:   ', num2str(Convergence(1, end))]);
    disp(['    Green: ', num2str(Convergence(2, end))]);
    disp(['    Blue:  ', num2str(Convergence(3, end))]);
    disp('The intensity: ');
    disp(['    Red:   ', num2str(intensity(1, :))]);
    disp(['    Green: ', num2str(intensity(2, :))]);
    disp(['    Blue:  ', num2str(intensity(3, :))]);
    disp(['The PSNR: ', num2str(psnr4mt)]);
    disp(['The SSIM: ', num2str(ssim4mt)]);
    disp(['The FSIM: ', num2str(fsim)]);
    disp(['The FSIMc: ', num2str(fsim_c)]);
end

%% Choose one result with the median indicators
choosen_In = [Indicator.(algorithm).SSIM];
median_index = find(choosen_In == median(choosen_In));
A = Indicator.(algorithm)(median_index(1));
A.Algorithm = algorithm;
A.PSNR_STD = std([Indicator.(algorithm).PSNR]);
A.SSIM_STD = std([Indicator.(algorithm).SSIM]);
A.FSIM_STD = std([Indicator.(algorithm).FSIM]);
A.Fitness_STD = std([Indicator.(algorithm).Fitness]')';
choose(i_algorithm) = A;
clear A;
figure();
mesh3('channel', Iout(:, :, 1), 'theme', 'red');
title('The red channel of the segmentation');
figure();
mesh3('channel', Iout(:, :, 2), 'theme', 'green');
title('The green channel of the segmentation');
figure();
mesh3('channel', Iout(:, :, 3), 'theme', 'blue');
title('The blue channel of the segmentation');
end
n = strsplit(img, '.');
save(cell2mat(strcat('MT_Indicator', '_', n(1), '_', num2str(level))), 'Indicator');
save(cell2mat(strcat('MT_Choose', '_', n(1), '_', num2str(level))), 'choose');
end
end
