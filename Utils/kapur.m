function fit = kapur(N, level, x, prob)
    fit = zeros(1, N);
    for i = 1:N
        PI0 = prob(1:x(i, 1));
        ind = PI0 == 0;
        ind = ind .* eps;
        PI0 = PI0 + ind;
        clear ind;
        w0 = sum(PI0);
        H0 = -sum((PI0 / w0) .* (log2(PI0 / w0)));
        fit(i) = fit(i) + H0;

        for j = 2:level - 1
            PI0 = prob(x(i, j - 1) + 1:x(i, j));
            ind = PI0 == 0;
            ind = ind .* eps;
            PI0 = PI0 + ind;
            clear ind;
            w0 = sum(PI0);
            H0 = -sum((PI0 / w0) .* (log2(PI0 / w0)));
            fit(i) = fit(i) + H0;
        end
    end
    fit = fit';
end
