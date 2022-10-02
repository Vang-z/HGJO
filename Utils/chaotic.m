function X = chaotic(N, dim, lb, ub)
    X = zeros(N, dim);
    w = 0.7;
    z = unifrnd(-1, 1, 1, dim);
    for i = 1:N
        z = sin((w * pi) ./ z);
        X(i, :) = lb + (1 + z) / 2 .* (ub - lb);
    end
    X = fix(X);
end
