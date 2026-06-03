%% TP3: Poroelasticidad - Problema de Terzaghi
clear; clc; close all

addpath('src/elements', 'src/assembly', 'src/solvers', 'src/postprocess', 'src/terzaghi_funs', 'meshes')

%% ---- Geometry ----
W  = 1.0;    % column width  [m]
L  = 10.0;   % column height [m]
nx = 2;      % elements in x (1 = pure 1D column)
ny = 20;     % elements in y (vertical refinement)
eleType = 'Q4';
npg     = 2;  % Gauss points per direction

mesh = generateColumnMesh(W, L, nx, ny, eleType);
nNod  = size(mesh.nodes, 1);
nElem = size(mesh.elements, 1);

fprintf('Mesh: %d nodes, %d elements (%s)\n', nNod, nElem, eleType);

%% ---- Applied load ----
sigma0 = 10000;   % [Pa]  uniform compressive traction on top face

%% =====================================================================
%% ---- TABLE 1 parameters ----
%% =====================================================================
params1.lambda = 4.0e7;       % [Pa]
params1.mu     = 4.0e7;       % [Pa]
params1.alpha  = 0.4;         % [-]
params1.Mbiot  = 4.1667e7;    % [Pa]
params1.phi    = 0.375;       % [-]
params1.kperm  = 1.01937e-9;  % [m^2]
params1.muf    = 1.0;         % [Pa*s]
params1.sigma0 = sigma0;
params1.L      = L;

Moed1  = params1.lambda + 2*params1.mu;
Kdr1   = Moed1;   % 1D oedometric = drained modulus
ratio1 = params1.alpha^2 * params1.Mbiot / Kdr1;
cv1    = (params1.kperm/params1.muf) / (1/params1.Mbiot + params1.alpha^2/Moed1);
p0_1   = params1.alpha * params1.Mbiot * sigma0 / (Moed1 + params1.alpha^2 * params1.Mbiot);

fprintf('\n=== Table 1 ===\n')
fprintf('  Moed        = %.4g Pa\n', Moed1)
fprintf('  alpha^2*M/Kdr = %.4f  (<<1 => converges)\n', ratio1)
fprintf('  cv          = %.4g m^2/s\n', cv1)
fprintf('  p0 (undrained) = %.4g Pa\n', p0_1)

%% ---- Build global matrices (Table 1) ----
[K1, H1, Q1, S1] = buildGlobalMatrices(mesh, eleType, params1, npg);
F1 = assembleLoad(mesh, eleType, sigma0, npg);

%% ---- Boundary conditions ----
% Mechanical: ux=0 on left (x=0) and right (x=W), uy=0 on bottom (y=0)
% Pressure:   p=0 on top (y=L)
[freeDofU1, freeDofP1] = buildCaseBC(mesh, W, L);

%% ---- Initial conditions (undrained) ----
% p(x,0) = p0 (uniform), u(x,0) from equilibrium
[U0_1, P0_1] = undrainedIC(K1, Q1, F1, freeDofU1, freeDofP1, p0_1);

%% ---- Task 1: Staggered, dt=1s, tmax=6000s ----
fprintf('\n--- Task 1: staggered, dt=1s ---\n')
dt1    = 1.0;
tmax1  = 6000.0;
tol    = 1e-8;
maxIter = 50;
saveTimes1 = [60, 300, 600, 1200, 3000, 6000];

[U1, P1, t1, iter1, ~] = staggeredSolver(...
    K1, H1, Q1, S1, ...
    F1, ...
    U0_1, P0_1, ...
    freeDofU1, freeDofP1, ...
    dt1, tmax1, tol, maxIter, saveTimes1);

fprintf('  Mean iterations/step: %.2f\n', mean(iter1(iter1>0)))

%% =====================================================================
%% ---- TABLE 2 parameters ----
%% =====================================================================
params2.lambda = 4.0e7;
params2.mu     = 4.0e7;
params2.alpha  = 0.4;
params2.Mbiot  = 6.06e9;
params2.phi    = 0.375;
params2.kperm  = 1.01937e-9;
params2.muf    = 1.0;
params2.sigma0 = sigma0;
params2.L      = L;

Moed2  = params2.lambda + 2*params2.mu;
Kdr2   = Moed2;
ratio2 = params2.alpha^2 * params2.Mbiot / Kdr2;
cv2    = (params2.kperm/params2.muf) / (1/params2.Mbiot + params2.alpha^2/Moed2);
p0_2   = params2.alpha * params2.Mbiot * sigma0 / (Moed2 + params2.alpha^2 * params2.Mbiot);

fprintf('\n=== Table 2 ===\n')
fprintf('  Moed        = %.4g Pa\n', Moed2)
fprintf('  alpha^2*M/Kdr = %.4f  (>>1 => diverges)\n', ratio2)
fprintf('  cv          = %.4g m^2/s\n', cv2)
fprintf('  p0 (undrained) = %.4g Pa\n', p0_2)

%% ---- Build global matrices (Table 2) ----
[K2, H2, Q2, S2] = buildGlobalMatrices(mesh, eleType, params2, npg);
F2 = assembleLoad(mesh, eleType, sigma0, npg);
[freeDofU2, freeDofP2] = buildCaseBC(mesh, W, L);
[U0_2, P0_2] = undrainedIC(K2, Q2, F2, freeDofU2, freeDofP2, p0_2);

%% ---- Task 2: Staggered Table 2, dt=0.1s, tmax=600s ----
fprintf('\n--- Task 2: staggered Table 2, dt=0.1s ---\n')
dt2    = 0.1;
tmax2  = 600.0;
saveTimes2 = [1, 5, 10, 30, 60, 300, 600];

[U2, P2, t2, iter2, pNorm2] = staggeredSolver(...
    K2, H2, Q2, S2, F2, U0_2, P0_2, freeDofU2, freeDofP2, ...
    dt2, tmax2, tol, maxIter, saveTimes2);

%% ---- Task 3: Fixed-stress split, Table 2 ----
fprintf('\n--- Task 3: fixed-stress split, Table 2, dt=0.1s ---\n')
% Contraction factor rho = alpha^2*M/(Kdr + alpha^2*M) ~ 0.89 => ~160 iters for 1e-8
% Use tol=1e-6 => ~120 iterations suffice
Smass2  = params2.Mbiot * S2;                        % int N'N dOmega
Lstab2  = (params2.alpha^2 / Kdr2) * Smass2;         % stabilization matrix
tol3    = 1e-6;
maxIter3 = 200;

[U3, P3, t3, iter3] = fixedStressSolver(...
    K2, H2, Q2, S2, F2, U0_2, P0_2, freeDofU2, freeDofP2, ...
    dt2, tmax2, tol3, maxIter3, saveTimes2, Lstab2);

fprintf('  Mean iterations/step: %.2f\n', mean(iter3(iter3>0)))

%% =====================================================================
%% ---- Report figures (saved to figs/) ----
%% =====================================================================
outdir = 'figs';
if ~exist(outdir, 'dir'); mkdir(outdir); end
set(groot, 'defaultAxesFontSize', 11, 'defaultLineLineWidth', 1.4, ...
           'defaultFigureColor', 'w');

% Pressure DOFs live on corner nodes (for Q4 = all nodes); take the x=0 column
[cn, nodeToP] = pressureNodeMap(mesh);
col0_n = cn(abs(mesh.nodes(cn,1)) < 1e-9);
[yP, iyp] = sort(mesh.nodes(col0_n,2));  col0_n = col0_n(iyp);
col0_p = nodeToP(col0_n);                          % rows in P_hist
allNod = (1:nNod)';
col0_u = allNod(abs(mesh.nodes(:,1)) < 1e-9);
[yU, iyu] = sort(mesh.nodes(col0_u,2));  col0_u = col0_u(iyu);
xq = linspace(0, L, 300)';

% --- Fig 1: Set 1 staggered — pressure + consolidation settlement increment ---
f = figure('Position',[100 100 760 360]);
cols = lines(numel(t1));
subplot(1,2,1); hold on; box on; grid on
for s = 1:numel(t1)
    pa = analyticalTerzaghi(xq, t1(s), params1);
    plot(pa/p0_1, xq, '--', 'Color', cols(s,:));
    plot(P1(col0_p,s)/p0_1, yP, 'o', 'Color', cols(s,:), 'MarkerSize', 4, ...
        'DisplayName', sprintf('%d s', t1(s)));
end
xlabel('p / p_0  [-]'); ylabel('Profundidad x  [m]'); title('(a) Presión'); ylim([0 L])
legend(findobj(gca,'Type','line','-and','Marker','o'),'Location','southeast','FontSize',8)
[~, ua0] = analyticalTerzaghi(xq, 0, params1);      % undrained reference
subplot(1,2,2); hold on; box on; grid on
for s = 1:numel(t1)
    [~, ua] = analyticalTerzaghi(xq, t1(s), params1);
    plot((ua-ua0)*1e6, xq, '--', 'Color', cols(s,:));
    plot((U1(2*col0_u,s)-U0_1(2*col0_u))*1e6, yU, 'o', 'Color', cols(s,:), 'MarkerSize', 4);
end
xlabel('u_y(t) - u_y(0)  [\mum]'); ylabel('Profundidad x  [m]');
title('(b) Asentamiento por consolidación'); ylim([0 L])
exportgraphics(f, fullfile(outdir,'tp3_tabla1.png'), 'Resolution', 200);

% --- Fig 2: Set 2 staggered — pressure-norm divergence ---
f = figure('Position',[100 100 430 330]); hold on; box on; grid on
for step = 1:5
    nrm = pNorm2{step};
    if numel(nrm) > 1
        semilogy(1:numel(nrm), nrm, '-o', 'MarkerSize', 4, ...
            'DisplayName', sprintf('paso %d (t=%.1f s)', step, step*dt2));
    end
end
set(gca,'YScale','log'); xlabel('Iteración k  [-]'); ylabel('||p^{(k)}||_2  [Pa]');
legend('Location','northwest','FontSize',8)
exportgraphics(f, fullfile(outdir,'tp3_divergencia.png'), 'Resolution', 200);

% --- Fig 3: Set 2 fixed-stress — pressure profiles ---
f = figure('Position',[100 100 430 360]); hold on; box on; grid on
cols = lines(numel(t3));
for s = 1:numel(t3)
    pa = analyticalTerzaghi(xq, t3(s), params2);
    plot(pa/p0_2, xq, '--', 'Color', cols(s,:));
    plot(P3(col0_p,s)/p0_2, yP, 'o', 'Color', cols(s,:), 'MarkerSize', 4, ...
        'DisplayName', sprintf('%g s', t3(s)));
end
xlabel('p / p_0  [-]'); ylabel('Profundidad x  [m]'); ylim([0 L])
legend(findobj(gca,'Type','line','-and','Marker','o'),'Location','southeast','FontSize',8)
exportgraphics(f, fullfile(outdir,'tp3_fixedstress.png'), 'Resolution', 200);

%% =====================================================================
%% ---- Error analysis & convergence studies (complement) ----
%% =====================================================================

% (1) Quantify the Task-1 match: max pressure error vs analytical (base mesh)
sIdx   = find(t1==300, 1);
pBase  = P1(col0_p, sIdx);
paBase = analyticalTerzaghi(yP, t1(sIdx), params1);
errBase = max(abs(pBase - paBase)) / p0_1 * 100;
fprintf('\n=== Error analysis & convergence ===\n')
fprintf('Base mesh %dx%d (ny x nx): max p-error vs analytical at t=%gs = %.3f %% of p0\n', ...
    ny, nx, t1(sIdx), errBase);

% (2) Delta-t study --------------------------------------------------------
% Set 1 (stable): accuracy vs dt on a fine mesh (100x10) -> backward Euler O(dt)
dtList = [2 1 0.5 0.25];
errDt  = zeros(size(dtList));
for i = 1:numel(dtList)
    errDt(i) = staggeredPError(W, L, 10, 100, eleType, npg, params1, p0_1, dtList(i), 300);
end
slopeDt = polyfit(log(dtList), log(errDt), 1);  slopeDt = slopeDt(1);
fprintf('\n-- Delta-t accuracy (Set 1, mesh 100x10, t=300s) --\n')
for i = 1:numel(dtList)
    fprintf('   dt = %5.3g s  ->  err = %.4f %% of p0\n', dtList(i), errDt(i));
end
fprintf('   temporal order (slope) ~ %.2f  (backward Euler -> 1)\n', slopeDt)

% Set 2 (unstable): divergence is dt-INDEPENDENT (rho = alpha^2 M/Kdr only)
fprintf('\n-- Delta-t effect on Set 2 staggered (diverges for all dt) --\n')
for dt = [0.2 0.1 0.05]
    ds = blowupStep(W, L, nx, ny, eleType, npg, params2, p0_2, dt, 5);
    fprintf('   dt = %5.3g s  ->  inner loop never converges; overflow at step %d (t=%.3g s)\n', ...
        dt, ds, ds*dt);
end

% (3) Mesh convergence -----------------------------------------------------
% To isolate the SPATIAL error we compare each mesh against a refined
% reference solution (400x40) computed with the SAME dt, so the temporal
% error is common to all meshes and cancels. (Comparing against the
% analytical solution instead lets the temporal floor mask the two finest
% meshes, hiding the spatial convergence.)
meshList  = [2 20; 10 100; 20 200];   % [nx ny]
dtMesh    = 0.5;  tEvalMesh = 30;     % dt cancels -> use a cheap, coarse dt
[yRef, pRef] = staggeredPProfile(W, L, 40, 400, eleType, npg, params1, dtMesh, tEvalMesh);
hMesh   = zeros(3,1);  errMesh = zeros(3,1);
for i = 1:3
    [yi, pi] = staggeredPProfile(W, L, meshList(i,1), meshList(i,2), eleType, npg, ...
                                 params1, dtMesh, tEvalMesh);
    errMesh(i) = max(abs(pi - interp1(yRef, pRef, yi))) / p0_1 * 100;
    hMesh(i)   = L / meshList(i,2);
end
slopeMesh = polyfit(log(hMesh), log(errMesh), 1);  slopeMesh = slopeMesh(1);
fprintf('\n-- Mesh convergence (Set 1, vs 400x40 reference, dt=%.2gs, t=%gs) --\n', dtMesh, tEvalMesh)
for i = 1:3
    fprintf('   mesh %3dx%-2d  h = %.3g m  ->  err = %.4f %% of p0\n', ...
        meshList(i,2), meshList(i,1), hMesh(i), errMesh(i));
end
fprintf('   spatial order (slope) ~ %.2f  (Q4 -> O(h^2))\n', slopeMesh)

% Convergence figure (temporal + spatial)
f = figure('Position',[100 100 760 320]);
subplot(1,2,1); loglog(dtList, errDt, 'o-', 'MarkerSize', 5); grid on; box on
xlabel('\Delta t  [s]'); ylabel('Error máx. de p  [% de p_0]');
title(sprintf('(a) Convergencia temporal  (pend. \\approx %.2f)', slopeDt));
subplot(1,2,2); loglog(hMesh, errMesh, 's-', 'MarkerSize', 6); grid on; box on
xlabel('h = L/n_y  [m]'); ylabel('Error máx. de p  [% de p_0]');
title(sprintf('(b) Convergencia espacial  (pend. \\approx %.2f)', slopeMesh));
exportgraphics(f, fullfile(outdir,'tp3_convergencia.png'), 'Resolution', 200);

% (4) Fixed-stress split convergence (Set 2) -------------------------------
% Same two studies for the stabilised scheme: temporal vs analytical on a
% fixed mesh, spatial vs a 400x40 reference at fixed dt.
dtListFs = [2 1 0.5 0.25];
errDtFs  = zeros(size(dtListFs));
for i = 1:numel(dtListFs)
    errDtFs(i) = fixedStressPError(W, L, 10, 100, eleType, npg, params2, p0_2, dtListFs(i), 30);
end
slopeDtFs = polyfit(log(dtListFs), log(errDtFs), 1);  slopeDtFs = slopeDtFs(1);
fprintf('\n-- Fixed-stress temporal (Set 2, mesh 100x10, t=30s) --\n')
for i = 1:numel(dtListFs)
    fprintf('   dt = %5.3g s  ->  err = %.4f %% of p0\n', dtListFs(i), errDtFs(i));
end
fprintf('   temporal order (slope) ~ %.2f\n', slopeDtFs)

[yRefFs, pRefFs] = fixedStressPProfile(W, L, 40, 400, eleType, npg, params2, 0.5, 30);
errMeshFs = zeros(3,1);
for i = 1:3
    [yi, pi] = fixedStressPProfile(W, L, meshList(i,1), meshList(i,2), eleType, npg, params2, 0.5, 30);
    errMeshFs(i) = max(abs(pi - interp1(yRefFs, pRefFs, yi))) / p0_2 * 100;
end
slopeMeshFs = polyfit(log(hMesh), log(errMeshFs), 1);  slopeMeshFs = slopeMeshFs(1);
fprintf('\n-- Fixed-stress spatial (Set 2, vs 400x40 reference, dt=0.5s, t=30s) --\n')
for i = 1:3
    fprintf('   mesh %3dx%-2d  h = %.3g m  ->  err = %.4f %% of p0\n', ...
        meshList(i,2), meshList(i,1), hMesh(i), errMeshFs(i));
end
fprintf('   spatial order (slope) ~ %.2f\n', slopeMeshFs)

f = figure('Position',[100 100 760 320]);
subplot(1,2,1); loglog(dtListFs, errDtFs, 'o-', 'MarkerSize', 5, 'Color', [0.83 0.33 0.10]);
grid on; box on
xlabel('\Delta t  [s]'); ylabel('Error máx. de p  [% de p_0]');
title(sprintf('(a) Convergencia temporal  (pend. \\approx %.2f)', slopeDtFs));
subplot(1,2,2); loglog(hMesh, errMeshFs, 's-', 'MarkerSize', 6, 'Color', [0.83 0.33 0.10]);
grid on; box on
xlabel('h = L/n_y  [m]'); ylabel('Error máx. de p  [% de p_0]');
title(sprintf('(b) Convergencia espacial  (pend. \\approx %.2f)', slopeMeshFs));
exportgraphics(f, fullfile(outdir,'tp3_convergencia_fs.png'), 'Resolution', 200);

% % --- Fig 4: problem schematic (compact; proportions not to scale) ---
% Ws = 1; Hs = 4;
% f = figure('Position',[100 100 280 320]); hold on; axis equal; axis off
% rectangle('Position',[0 0 Ws Hs],'FaceColor',[0.86 0.9 0.96],'EdgeColor','k','LineWidth',1.2);
% plot([0 Ws],[Hs+0.5 Hs+0.5],'k','LineWidth',1.2);
% xa = linspace(0.1, Ws-0.1, 7);
% for k = 1:numel(xa)
%     plot([xa(k) xa(k)],[Hs+0.5 Hs+0.12],'k','LineWidth',1.0);
%     plot(xa(k), Hs+0.12, 'kv','MarkerFaceColor','k','MarkerSize',4);
% end
% text(Ws/2, Hs+0.9, '\sigma_0 = 10000 Pa','HorizontalAlignment','center','FontSize',10);
% text(Ws/2, Hs-0.30, 'p = 0','HorizontalAlignment','center','FontSize',9,'Color',[0 0 0.7]);
% hs = Hs/8;
% for yy = 0:hs:Hs-hs, plot([-0.14 0.02],[yy yy+0.16],'k','LineWidth',0.4); end
% for yy = 0:hs:Hs-hs, plot([Ws-0.02 Ws+0.14],[yy yy+0.16],'k','LineWidth',0.4); end
% for xx = 0:Ws/8:Ws-Ws/8, plot([xx xx+Ws/8],[-0.14 0.02],'k','LineWidth',0.4); end
% plot([0 0],[0 Hs],'k','LineWidth',1.2); plot([Ws Ws],[0 Hs],'k','LineWidth',1.2);
% plot([0 Ws],[0 0],'k','LineWidth',1.2);
% text(-0.30, Hs/2, 'u_x = 0','FontSize',9,'Rotation',90,'HorizontalAlignment','center');
% text(Ws/2, -0.45, 'u_y = 0  (impermeable)','FontSize',9,'HorizontalAlignment','center');
% text(Ws+0.52, Hs/2, 'L = 10 m','FontSize',9,'Rotation',90,'HorizontalAlignment','center');
% plot([Ws+0.28 Ws+0.28],[0 Hs/3],'k','LineWidth',1); plot(Ws+0.28, Hs/3, 'k^','MarkerFaceColor','k','MarkerSize',4);
% text(Ws+0.36, Hs/3-0.12, 'x','FontSize',11);
% xlim([-0.7 2.0]); ylim([-0.8 Hs+1.4]);
% exportgraphics(f, fullfile(outdir,'tp3_esquema.png'), 'Resolution', 200);
% 
% fprintf('Report figures saved to %s/\n', outdir);
%
% fprintf('\nDone.\n')

%% ---- local helpers (convergence studies) ----
function eP = staggeredPError(W, L, nx, ny, eleType, npg, par, p0, dt, tEval)
% Max relative pressure error [% of p0] vs the analytical solution at tEval,
% for a staggered Set-1-type (stable) run on an nx-by-ny mesh.
    m = generateColumnMesh(W, L, nx, ny, eleType);
    [K, H, Q, S] = buildGlobalMatrices(m, eleType, par, npg);
    F = assembleLoad(m, eleType, par.sigma0, npg);
    [fU, fP] = buildCaseBC(m, W, L);
    [U0, P0] = undrainedIC(K, Q, F, fU, fP, p0);
    [~, Ph, th] = staggeredSolver(K, H, Q, S, F, U0, P0, fU, fP, dt, tEval, 1e-8, 50, tEval);
    [c, n2p] = pressureNodeMap(m);
    c0 = c(abs(m.nodes(c,1)) < 1e-9);          % corner nodes on x=0 column
    yv = m.nodes(c0, 2);
    pf = Ph(n2p(c0), end);
    pa = analyticalTerzaghi(yv, th(end), par);
    eP = max(abs(pf - pa)) / p0 * 100;
end

function [yv, pv] = staggeredPProfile(W, L, nx, ny, eleType, npg, par, dt, tEval)
% Centreline pressure profile (sorted by height) of a stable staggered run,
% used for mesh-to-reference spatial convergence (temporal error cancels).
    m = generateColumnMesh(W, L, nx, ny, eleType);
    [K, H, Q, S] = buildGlobalMatrices(m, eleType, par, npg);
    F = assembleLoad(m, eleType, par.sigma0, npg);
    [fU, fP] = buildCaseBC(m, W, L);
    p0 = par.alpha*par.Mbiot*par.sigma0 / ((par.lambda+2*par.mu) + par.alpha^2*par.Mbiot);
    [U0, P0] = undrainedIC(K, Q, F, fU, fP, p0);
    [~, Ph] = staggeredSolver(K, H, Q, S, F, U0, P0, fU, fP, dt, tEval, 1e-8, 50, tEval);
    [c, n2p] = pressureNodeMap(m);
    c0 = c(abs(m.nodes(c,1)) < 1e-9);
    [yv, o] = sort(m.nodes(c0, 2));
    pv = Ph(n2p(c0), end);  pv = pv(o);
end

function eP = fixedStressPError(W, L, nx, ny, eleType, npg, par, p0, dt, tEval)
% Max relative pressure error [% of p0] vs analytical for a fixed-stress run.
    m = generateColumnMesh(W, L, nx, ny, eleType);
    [K, H, Q, S] = buildGlobalMatrices(m, eleType, par, npg);
    F = assembleLoad(m, eleType, par.sigma0, npg);
    [fU, fP] = buildCaseBC(m, W, L);
    [U0, P0] = undrainedIC(K, Q, F, fU, fP, p0);
    Lstab = (par.alpha^2 / (par.lambda+2*par.mu)) * par.Mbiot * S;
    [~, Ph, th] = fixedStressSolver(K, H, Q, S, F, U0, P0, fU, fP, dt, tEval, 1e-6, 300, tEval, Lstab);
    [c, n2p] = pressureNodeMap(m);
    c0 = c(abs(m.nodes(c,1)) < 1e-9);
    yv = m.nodes(c0, 2);
    pf = Ph(n2p(c0), end);
    pa = analyticalTerzaghi(yv, th(end), par);
    eP = max(abs(pf - pa)) / p0 * 100;
end

function [yv, pv] = fixedStressPProfile(W, L, nx, ny, eleType, npg, par, dt, tEval)
% Centreline pressure profile of a fixed-stress run, for the spatial study.
    m = generateColumnMesh(W, L, nx, ny, eleType);
    [K, H, Q, S] = buildGlobalMatrices(m, eleType, par, npg);
    F = assembleLoad(m, eleType, par.sigma0, npg);
    [fU, fP] = buildCaseBC(m, W, L);
    p0 = par.alpha*par.Mbiot*par.sigma0 / ((par.lambda+2*par.mu) + par.alpha^2*par.Mbiot);
    [U0, P0] = undrainedIC(K, Q, F, fU, fP, p0);
    Lstab = (par.alpha^2 / (par.lambda+2*par.mu)) * par.Mbiot * S;
    [~, Ph] = fixedStressSolver(K, H, Q, S, F, U0, P0, fU, fP, dt, tEval, 1e-6, 300, tEval, Lstab);
    [c, n2p] = pressureNodeMap(m);
    c0 = c(abs(m.nodes(c,1)) < 1e-9);
    [yv, o] = sort(m.nodes(c0, 2));
    pv = Ph(n2p(c0), end);  pv = pv(o);
end

function ds = blowupStep(W, L, nx, ny, eleType, npg, par, p0, dt, tmax)
% First time step at which the staggered solution overflows (non-finite),
% used to show the dt-independence of the Set-2 divergence.
    m = generateColumnMesh(W, L, nx, ny, eleType);
    [K, H, Q, S] = buildGlobalMatrices(m, eleType, par, npg);
    F = assembleLoad(m, eleType, par.sigma0, npg);
    [fU, fP] = buildCaseBC(m, W, L);
    [U0, P0] = undrainedIC(K, Q, F, fU, fP, p0);
    [~, ~, ~, ~, pN] = staggeredSolver(K, H, Q, S, F, U0, P0, fU, fP, dt, tmax, 1e-8, 50, tmax);
    ds = numel(pN);
    for s = 1:numel(pN)
        v = pN{s};
        if ~isempty(v) && any(~isfinite(v)); ds = s; break; end
    end
end
