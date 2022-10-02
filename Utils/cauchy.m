function p = cauchy(n, x0, gamma)
    p =	gamma ./ (pi * (gamma.^2 + (rand(n) - x0).^2));
end
