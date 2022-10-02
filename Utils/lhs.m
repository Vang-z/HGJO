function X = lhs(N, dim, lb, ub)
    % Use Latin Hypercube Sampling to initialize the populations
    X = lb + lhsdesign(N, dim) .* (ub - lb);
end
