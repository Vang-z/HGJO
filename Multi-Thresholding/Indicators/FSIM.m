function [FSIM, FSIMc] = FSIM(Iin, Iout)
%% FSIM: Calculate the Feature Similarity index metric to check the similarity between two images based on their internal feather.
%% Dimension: n_var --- dimensions of decision space

%% Input:
%                      Dimension                                      Description
%      Iin             row x col                                      The raw Image
%      Iout            row x col                                      The output Image

%% Output:
%                      Dimension                                      Description
%      FSIM            1 x 1                                          FSIM compares two images based on their luminance components only
%      FSIMc           1 x 1                                          FSIMc considers the chromatic information in addition to the luminance

%% Reference
%
%      Title:   FSIM: A Feature Similarity Index for Image Quality Assessment.
%      Authors: Lin Zhang; Lei Zhang; Xuanqin Mou; David Zhang.
%      DOI:     https://doi.org/10.1109/TIP.2011.2109730

%% NOTICE:
%
%  I acknowledge that this code has been written using a large portion of the following Paper:
%
%      Title:   FSIM: A Feature Similarity Index for Image Quality Assessment.
%      Authors: Lin Zhang; Lei Zhang; Xuanqin Mou; David Zhang.
%      DOI:     https://doi.org/10.1109/TIP.2011.2109730

[row, col] = size(Iin(:,:,1));
I1 = ones(row, col);
I2 = ones(row, col);
Q1 = ones(row, col);
Q2 = ones(row, col);

if size(Iin, 3) == 3
    Y1 = 0.299 * double(Iin(:, :, 1)) + 0.587 * double(Iin(:, :, 2)) + 0.114 * double(Iin(:, :, 3));
    Y2 = 0.299 * double(Iout(:, :, 1)) + 0.587 * double(Iout(:, :, 2)) + 0.114 * double(Iout(:, :, 3));
    I1 = 0.596 * double(Iin(:, :, 1)) - 0.274 * double(Iin(:, :, 2)) - 0.322 * double(Iin(:, :, 3));
    I2 = 0.596 * double(Iout(:, :, 1)) - 0.274 * double(Iout(:, :, 2)) - 0.322 * double(Iout(:, :, 3));
    Q1 = 0.211 * double(Iin(:, :, 1)) - 0.523 * double(Iin(:, :, 2)) + 0.312 * double(Iin(:, :, 3));
    Q2 = 0.211 * double(Iout(:, :, 1)) - 0.523 * double(Iout(:, :, 2)) + 0.312 * double(Iout(:, :, 3));
else
    Y1 = Iin;
    Y2 = Iout;
end

Y1 = double(Y1);
Y2 = double(Y2);
minDimension = min(row, col);
F = max(1,round(minDimension / 256));
aveKernel = fspecial('average', F);

aveI1 = conv2(I1, aveKernel, 'same');
aveI2 = conv2(I2, aveKernel, 'same');
I1 = aveI1(1:F:row, 1:F:col);
I2 = aveI2(1:F:row, 1:F:col);

aveQ1 = conv2(Q1, aveKernel, 'same');
aveQ2 = conv2(Q2, aveKernel, 'same');
Q1 = aveQ1(1:F:row, 1:F:col);
Q2 = aveQ2(1:F:row, 1:F:col);

aveY1 = conv2(Y1, aveKernel, 'same');
aveY2 = conv2(Y2, aveKernel, 'same');
Y1 = aveY1(1:F:row, 1:F:col);
Y2 = aveY2(1:F:row, 1:F:col);

PC1 = phasecong2(Y1);
PC2 = phasecong2(Y2);

dx = [3, 0, -3; 10, 0, -10; 3, 0, -3] / 16;
dy = [3, 10, 3; 0, 0, 0; -3, -10, -3] / 16;
IxY1 = conv2(Y1, dx, 'same');
IyY1 = conv2(Y1, dy, 'same');
gradientMap1 = sqrt(IxY1.^2 + IyY1.^2);

IxY2 = conv2(Y2, dx, 'same');
IyY2 = conv2(Y2, dy, 'same');
gradientMap2 = sqrt(IxY2.^2 + IyY2.^2);

T1 = 0.85;
T2 = 160;
PCSimMatrix = (2 * PC1 .* PC2 + T1) ./ (PC1.^2 + PC2.^2 + T1);
gradientSimMatrix = (2 * gradientMap1.*gradientMap2 + T2) ./ (gradientMap1.^2 + gradientMap2.^2 + T2);
PCm = max(PC1, PC2);
SimMatrix = gradientSimMatrix .* PCSimMatrix .* PCm;
FSIM = sum(sum(SimMatrix)) / sum(sum(PCm));

T3 = 200;
T4 = 200;
ISimMatrix = (2 * I1 .* I2 + T3) ./ (I1.^2 + I2.^2 + T3);
QSimMatrix = (2 * Q1 .* Q2 + T4) ./ (Q1.^2 + Q2.^2 + T4);

lambda = 0.03;

SimMatrixC = gradientSimMatrix .* PCSimMatrix .* real((ISimMatrix .* QSimMatrix).^lambda) .* PCm;
FSIMc = sum(sum(SimMatrixC)) / sum(sum(PCm));
end


function [ResultPC] = phasecong2(im)
    % Number of wavelet scales.
    nscale = 4;
    % Number of filter orientations.
    norient = 4;
    % Wavelength of smallest scale filter.
    minWaveLength = 6;
    % Scaling factor between successive filters.
    mult = 2;
    % Ratio of the standard deviation of the Gaussian describing the log Gabor filter's transfer function in
    % the frequency domain to the filter center frequency.
    sigmaOnf = 0.55;
    % Ratio of angular interval between filter orientations and the standard deviation of the angular Gaussian
    % function used to construct filters in the freq. plane.
    dThetaOnSigma = 1.2;
    % No of standard deviations of the noise energy beyond the mean at which we set the noise threshold point.
    % below which phase congruency values get penalized.
    k = 2.0;
    % Used to prevent division by zero.
    epsilon = 0.0001;
    % Calculate the standard deviation of the angular Gaussian function used to construct filters in the freq. plane.
    thetaSigma = pi / norient / dThetaOnSigma;

    [row,col] = size(im);
    % Fourier transform of image
    imagefft = fft2(im);

    zero = zeros(row, col);
    % Array of convolution results.
    EO = cell(nscale, norient);

    estMeanE2n = [];
    % Array of inverse FFTs of filters
    ifftFilterArray = cell(1, nscale);

    % Pre-compute some stuff to speed up filter construction
    % Set up X and Y matrices with ranges normalised to +/- 0.5
    % The following code adjusts things appropriately for odd and even values of row and columns.
    if mod(col, 2)
        xrange = [-(col - 1) / 2:(col - 1) / 2] / (col - 1);
    else
        xrange = [-col / 2:(col / 2 - 1)] / col;
    end

    if mod(row, 2)
        yrange = [-(row - 1) / 2:(row - 1) / 2] / (row - 1);
    else
        yrange = [-row / 2:(row / 2 - 1)] / row;
    end

    [x,y] = meshgrid(xrange, yrange);
    % Matrix values contain normalised radius from centre.
    radius = sqrt(x.^2 + y.^2);
    % Matrix values contain polar angle. (note -ve y is used to give +ve anti-clockwise angles)
    theta = atan2(-y, x);

    % Quadrant shift radius and theta so that filters
    radius = ifftshift(radius);
    % are constructed with 0 frequency at the corners.
    theta = ifftshift(theta);
    % Get rid of the 0 radius value at the 0 frequency point (now at top-left corner)
    % so that taking the log of the radius will not cause trouble.
    radius(1,1) = 1;

    sintheta = sin(theta);
    costheta = cos(theta);
    clear x; clear y; clear theta;

    % Filters are constructed in terms of two components.
    % 1) The radial component, which controls the frequency band that the filter
    %    responds to
    % 2) The angular component, which controls the orientation that the filter
    %    responds to.
    % The two components are multiplied together to construct the overall filter.

    % Construct the radial filter components...

    % First construct a low-pass filter that is as large as possible, yet falls
    % away to zero at the boundaries.  All log Gabor filters are multiplied by
    % this to ensure no extra frequencies at the 'corners' of the FFT are
    % incorporated as this seems to upset the normalisation process when
    % calculating phase congrunecy.
    % Radius .45, 'sharpness' 15
    lp = lowpassfilter([row, col], 0.45, 15);

    logGabor = cell(1, nscale);

    for s = 1:nscale
        wavelength = minWaveLength*mult^(s - 1);
        % Centre frequency of filter.
        fo = 1.0 / wavelength;
        logGabor{s} = exp((-(log(radius / fo)).^2) / (2 * log(sigmaOnf)^2));
        % Apply low-pass filter
        logGabor{s} = logGabor{s} .* lp;
        % Set the value at the 0 frequency point of the filter back to zero (undo the radius fudge).
        logGabor{s}(1, 1) = 0;
    end

    % Then construct the angular filter components...
    spread = cell(1, norient);

    for o = 1:norient
        % Filter angle.
        angl = (o - 1) * pi / norient;
        % For each point in the filter matrix calculate the angular distance from
        % the specified filter orientation.  To overcome the angular wrap-around
        % problem sine difference and cosine difference values are first computed
        % and then the atan2 function is used to determine angular distance.
        % Difference in sine.
        ds = sintheta * cos(angl) - costheta * sin(angl);
        % Difference in cosine.
        dc = costheta * cos(angl) + sintheta * sin(angl);
        % Absolute angular distance.
        dtheta = abs(atan2(ds, dc));
        % Calculate the angular filter component.
        spread{o} = exp((-dtheta.^2) / (2 * thetaSigma^2));
    end

    % The main loop...
    EnergyAll(row,col) = 0;
    AnAll(row,col) = 0;

    for o = 1:norient
        % Initialize accumulator matrices.
        sumE_ThisOrient = zero;
        sumO_ThisOrient = zero;
        sumAn_ThisOrient = zero;
        Energy = zero;
        for s = 1:nscale
            % Multiply radial and angular components to get the filter.
            filter = logGabor{s} .* spread{o};
            % Note rescaling to match power
            ifftFilt = real(ifft2(filter)) * sqrt(row * col);
            % record ifft2 of filter
            ifftFilterArray{s} = ifftFilt;
            % Convolve image with even and odd filters returning the result in EO
            EO{s, o} = ifft2(imagefft .* filter);

            % Amplitude of even & odd filter response.
            An = abs(EO{s, o});
            % Sum of amplitude responses.
            sumAn_ThisOrient = sumAn_ThisOrient + An;
            % Sum of even filter convolution results.
            sumE_ThisOrient = sumE_ThisOrient + real(EO{s, o});
            % Sum of odd filter convolution results.
            sumO_ThisOrient = sumO_ThisOrient + imag(EO{s, o});
            % Record mean squared filter value at smallest
            if s == 1
                % scale. This is used for noise estimation.
                EM_n = sum(sum(filter.^2));
                % Record the maximum An over all scales.
                maxAn = An;
            else
                maxAn = max(maxAn, An);
            end
        end

        % Get weighted mean filter response vector, this gives the weighted mean phase angle.
        XEnergy = sqrt(sumE_ThisOrient.^2 + sumO_ThisOrient.^2) + epsilon;
        MeanE = sumE_ThisOrient ./ XEnergy;
        MeanO = sumO_ThisOrient ./ XEnergy;

        % Now calculate An(cos(phase_deviation) - | sin(phase_deviation)) | by using dot and cross products between
        % the weighted mean filter response vector and the individual filter response vectors at each scale.  This
        % quantity is phase congruency multiplied by An, which we call energy.
        for s = 1:nscale
          E = real(EO{s, o});
          O = imag(EO{s, o});
          Energy = Energy + E .* MeanE + O .* MeanO - abs(E .* MeanO - O .* MeanE);
        end

        % Compensate for noise We estimate the noise power from the energy squared response at the smallest scale.
        % If the noise is Gaussian the energy squared will have a Chi-squared 2DOF pdf.  We calculate the median
        % energy squared response  as this is a robust statistic. From this we estimate the mean. The estimate of
        % noise power is obtained by dividing the mean squared energy value by the mean squared filter value.

        medianE2n = median(reshape(abs(EO{1, o}).^2, 1, row * col));
        meanE2n = -medianE2n / log(0.5);
        estMeanE2n(o) = meanE2n;
        % Estimate of noise power.
        noisePower = meanE2n / EM_n;

        % Now estimate the total energy^2 due to noise Estimate for sum(An^2) + sum(Ai.*Aj.*(cphi.*cphj + sphi.*sphj))
        EstSumAn2 = zero;
        for s = 1:nscale
            EstSumAn2 = EstSumAn2 + ifftFilterArray{s}.^2;
        end

        EstSumAiAj = zero;
        for si = 1:(nscale - 1)
            for sj = (si + 1):nscale
                EstSumAiAj = EstSumAiAj + ifftFilterArray{si} .* ifftFilterArray{sj};
            end
        end
        sumEstSumAn2 = sum(sum(EstSumAn2));
        sumEstSumAiAj = sum(sum(EstSumAiAj));

        EstNoiseEnergy2 = 2 * noisePower*sumEstSumAn2 + 4 * noisePower*sumEstSumAiAj;
        % Rayleigh parameter
        tau = sqrt(EstNoiseEnergy2 / 2);
        % Expected value of noise energy
        EstNoiseEnergy = tau * sqrt(pi / 2);
        EstNoiseEnergySigma = sqrt((2 - pi / 2) * tau^2);
        % Noise threshold
        T =  EstNoiseEnergy + k * EstNoiseEnergySigma;

        % The estimated noise effect calculated above is only valid for the PC_1 measure. The PC_2 measure does not
        % lend itself readily to the same analysis. However empirically it seems that the noise effect is overestimated
        % roughly by a factor of 1.7 for the filter parameters used here.

        % Empirical rescaling of the estimated noise effect to uit the PC_2 phase congruency measure
        T = T / 1.7;
        % Apply noise threshold
        Energy = max(Energy - T, zero);
        EnergyAll = EnergyAll + Energy;
        AnAll = AnAll + sumAn_ThisOrient;
    end
    ResultPC = EnergyAll ./ AnAll;
end


function f = lowpassfilter(sze, cutoff, n)
    if cutoff < 0 || cutoff > 0.5
	    error('cutoff frequency must be between 0 and 0.5');
    end

    if rem(n,1) ~= 0 || n < 1
	    error('n must be an integer >= 1');
    end

    if length(sze) == 1
	    row = sze;
	    col = sze;
    else
	    row = sze(1);
	    col = sze(2);
    end

    % Set up X and Y matrices with ranges normalised to +/- 0.5
    % The following code adjusts things appropriately for odd and even values of row and columns.
    if mod(col, 2)
	    xrange = [-(col - 1) / 2:(col - 1) / 2] / (col - 1);
    else
	    xrange = [-col / 2:(col / 2 - 1)] / col;
    end

    if mod(row, 2)
	    yrange = [-(row - 1) / 2:(row - 1) / 2] / (row - 1);
    else
	    yrange = [-row / 2:(row / 2 - 1)] / row;
    end

    [x, y] = meshgrid(xrange, yrange);
    % A matrix with every pixel = radius relative to centre.
    radius = sqrt(x.^2 + y.^2);
    f = ifftshift(1 ./ (1.0 + (radius ./ cutoff).^(2 * n)));
end
