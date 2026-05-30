function plotPressureNorm(pNormHist, dt, titleStr)
    % Find first step with more than 1 iteration
    nSteps = numel(pNormHist);
    figure('Name', titleStr);
    hold on
    for step = 1:min(nSteps, 5)
        nrm = pNormHist{step};
        if numel(nrm) > 1
            plot(1:numel(nrm), nrm, '-o', 'DisplayName', sprintf('step %d (t=%.2gs)', step, step*dt))
        end
    end
    xlabel('Iteration k');
    ylabel('||p^{(k)}||_2  [Pa]');
    title(titleStr);
    legend show
    grid on
    hold off
end
