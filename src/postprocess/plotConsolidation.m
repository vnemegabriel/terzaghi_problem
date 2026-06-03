function plotConsolidation(mesh, U_hist, P_hist, t_hist, params, titleStr)
% plotConsolidation  Compare FEM results vs analytical solution.
%   Extracts vertical displacement at x=0.5m and pressure at x=0.5m
%   along the column height for each saved time.

L      = params.L;
sigma0 = params.sigma0;
nSave  = numel(t_hist);
colors = lines(nSave);

% Query points along column centreline (x_horiz=0.5m)
x_query = linspace(0, L, 200)';
Moed = params.lambda + 2*params.mu;
Mbiot = params.Mbiot;
alpha = params.alpha;
p0   = alpha * Mbiot * sigma0 / (Moed + alpha^2 * Mbiot);

% Pressure DOFs live on corner nodes only (mixed u-p formulation).
cornerNodes = pressureNodeMap(mesh);
yP = mesh.nodes(cornerNodes, 2);

% --- Pressure plot ---
figure('Name', ['Pressure - ' titleStr]);
hold on
for s = 1:nSave
    % FEM: extract pressure at nodes along y-axis (x~0 column)
    [p_anal, ~] = analyticalTerzaghi(x_query, t_hist(s), params);
    plot(p_anal / p0, x_query, '--', 'Color', colors(s,:), 'LineWidth', 1.2)

    % FEM nodal pressures at corner nodes, sorted by height
    pFEM_nodes = P_hist(:, s);
    [ySorted, idx] = sort(yP);
    pSorted = pFEM_nodes(idx);
    plot(pSorted / p0, ySorted, '-o', 'Color', colors(s,:), 'MarkerSize', 3, ...
        'DisplayName', sprintf('t=%.0f s', t_hist(s)))
end
xlabel('p / p_0  [-]');
ylabel('x [m]');
title(['Pore pressure: ' titleStr]);
legend show
grid on
hold off

% --- Displacement plot ---
figure('Name', ['Displacement - ' titleStr]);
hold on
for s = 1:nSave
    [~, u_anal] = analyticalTerzaghi(x_query, t_hist(s), params);
    plot(u_anal * 1e3, x_query, '--', 'Color', colors(s,:), 'LineWidth', 1.2)

    uFEM_all = reshape(U_hist(:,s), 2, [])';  % [nNod x 2]
    yNodes   = mesh.nodes(:,2);
    [ySorted, idx] = sort(yNodes);
    uySorted = uFEM_all(idx, 2);
    plot(uySorted * 1e3, ySorted, 'o', 'Color', colors(s,:), 'MarkerSize', 6, ...
        'DisplayName', sprintf('t=%.0f s', t_hist(s)))
end
xlabel('u_y [mm]');
ylabel('x [m]');
title(['Displacement: ' titleStr]);
legend show
grid on
hold off
end
