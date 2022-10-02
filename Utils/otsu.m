function fit = otsu(N, level, x, prob)
    fit = zeros(1, N);
    for j = 1:N
        fit(j) = sum(prob(1:x(j, 1) - 1)) * (sum((1:x(j, 1) - 1) .* prob(1:x(j, 1) - 1) / ...
                 sum(prob(1:x(j, 1) - 1))) - sum((1:255) .* prob(1:255)) ) ^ 2;
        for jlevel = 2:level - 1
            fit(j) = fit(j) + sum(prob(x(j, jlevel - 1):x(j, jlevel) - 1)) * ...
                     (sum((x(j, jlevel - 1):x(j, jlevel) - 1) .* prob(x(j, jlevel - 1):x(j, jlevel) - 1) / ...
                     sum(prob(x(j, jlevel - 1):x(j, jlevel) - 1))) - sum((1:255) .* prob(1:255))) ^ 2;
        end
        fit(j) = fit(j) + sum(prob(x(j, level - 1):255)) * (sum((x(j, level - 1):255) .* ...
                 prob(x(j, level - 1):255) / sum(prob(x(j, level - 1):255))) - ...
                 sum((1:255) .* prob(1:255))) ^ 2;
        if isnan(fit(j)) || isinf(fit(j))
            fit(j) = eps;
        end
    end
    fit = fit';
end
