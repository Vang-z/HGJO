data = Table(3:5:end, :);
[~, ~, rk] = friedman(data, 1);
rk = rk.meanranks
