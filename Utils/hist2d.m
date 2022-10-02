function hist2d(Iin)
    % Compute the "histogram"
    % Allocate a histogram of 256 gray levels by 256 gray levels.
%    figure();
%    [rows, columns] = size(Iin);
%    windowWidth = 21;
%    kernel = ones(windowWidth)/windowWidth^2;
%    Bin = imfilter(Iin, kernel);
%    hist = zeros(256, 256);
%    for row = 1:rows
%        for column = 1 : columns
%            index1 = Iin(row, column);
%            index2 = round(Bin(row, column));
%            hist(index1 + 1, index2 + 1) = hist(index1 + 1, index2 + 1) + 1;
%        end
%    end
%    meshc(hist, 'FaceColor', 'default', 'EdgeColor', 'none');
%    set(gca, 'xaxislocation', 'top', 'yaxislocation', 'left', 'ydir', 'reverse');
%    colormap('turbo');
%    colorbar;
%    hold on;
%    fmesh(0, 'LineStyle', '-.', 'EdgeColor', '#6D597A', 'FaceAlpha', 0);
%    grid on;
%    title('The 2D Histogram');

    figure();
    [n_countR, x_valueR] = imhist(Iin(:, :, 1));
    [n_countG, x_valueG] = imhist(Iin(:, :, 2));
    [n_countB, x_valueB] = imhist(Iin(:, :, 3));
    n_countR(1) = 0;
    n_countG(1) = 0;
    n_countB(1) = 0;
    n_countR(256) = 0;
    n_countG(256) = 0;
    n_countB(256) = 0;
    fill(x_valueR, n_countR, 'r');
    hold on;
    fill(x_valueG, n_countG, 'g');
    hold on;
    fill(x_valueB, n_countB, 'b');
    hold on;
    title('The grayscale map');
end
