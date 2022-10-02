function plotf(lb, ub, f_obj)
    x_label = lb:2:ub;
    y_label = x_label;
    L = length(x_label);
    for i = 1:L
        for j = 1:L
            f(i, j) = f_obj([x_label(i), y_label(j)]);
        end
    end
    surfc(x_label, y_label, f);
    colormap('turbo');
    grid on;
end
