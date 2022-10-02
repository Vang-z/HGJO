function plot4Res(choose)
    c = {choose.Convergence};
    % plot Convergence
    figure;
    for i = 1:length(c)
        Convergence = cell2mat(c(i));
        semilogy(Convergence(1, :), 'linewidth', 1.5);
        hold on;
    end
    title('Objective space')
    xlabel('Iteration');
    ylabel('Best score obtained so far');
    axis tight;
    grid on;
    box on;
    legend({choose.Algorithm});
end
