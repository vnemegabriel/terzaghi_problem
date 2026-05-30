%% TP3: Poroelasticidad - Problema de Terzaghi
clear; clc; close all

addpath('src/elements', 'src/assembly', 'src/solvers', 'src/postprocess', 'src/terzaghi_funs', 'src/utils', 'meshes')

%% ---- Geometry ----
W  = 1.0;    % column width  [m]
L  = 10.0;   % column height [m]
nx = 2;      % elements in x (1 = pure 1D column)
ny = 20;     % elements in y (vertical refinement)
eleType = 'Q8';
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
p1.lambda = 4.0e7;       % [Pa]
p1.mu     = 4.0e7;       % [Pa]
p1.alpha  = 0.4;         % [-]
p1.Mbiot  = 4.1667e7;    % [Pa]
p1.phi    = 0.375;       % [-]
p1.kperm  = 1.01937e-9;  % [m^2]
p1.muf    = 1.0;         % [Pa*s]
p1.sigma0 = sigma0;
p1.L      = L;

Moed1  = p1.lambda + 2*p1.mu;
Kdr1   = Moed1;   % 1D oedometric = drained modulus
ratio1 = p1.alpha^2 * p1.Mbiot / Kdr1;
cv1    = (p1.kperm/p1.muf) / (1/p1.Mbiot + p1.alpha^2/Moed1);
p0_1   = p1.alpha * p1.Mbiot * sigma0 / (Moed1 + p1.alpha^2 * p1.Mbiot);

fprintf('\n=== Table 1 ===\n')
fprintf('  Moed        = %.4g Pa\n', Moed1)
fprintf('  alpha^2*M/Kdr = %.4f  (<<1 => converges)\n', ratio1)
fprintf('  cv          = %.4g m^2/s\n', cv1)
fprintf('  p0 (undrained) = %.4g Pa\n', p0_1)

%% ---- Build global matrices (Table 1) ----
[K1, H1, Q1, S1] = buildGlobalMatrices(mesh, eleType, p1, npg);
F1 = assembleLoad(mesh, eleType, sigma0, npg);

%% ---- Boundary conditions ----
% Mechanical: ux=0 on left (x=0) and right (x=W), uy=0 on bottom (y=0)
% Pressure:   p=0 on top (y=L)
[freeDofU1, freeDofP1] = buildCaseBC(mesh, W, L);

%% ---- Initial conditions (undrained) ----
% p(x,0) = p0 (uniform), u(x,0) from equilibrium
[U0_1, P0_1] = undrainedIC(K1, Q1, F1, freeDofU1, freeDofP1, p0_1, nNod);

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
plotConsolidation(mesh, U1, P1, t1, p1, 'Table 1 - Staggered')

%% =====================================================================
%% ---- TABLE 2 parameters ----
%% =====================================================================
p2.lambda = 4.0e7;
p2.mu     = 4.0e7;
p2.alpha  = 0.4;
p2.Mbiot  = 6.06e9;
p2.phi    = 0.375;
p2.kperm  = 1.01937e-9;
p2.muf    = 1.0;
p2.sigma0 = sigma0;
p2.L      = L;

Moed2  = p2.lambda + 2*p2.mu;
Kdr2   = Moed2;
ratio2 = p2.alpha^2 * p2.Mbiot / Kdr2;
cv2    = (p2.kperm/p2.muf) / (1/p2.Mbiot + p2.alpha^2/Moed2);
p0_2   = p2.alpha * p2.Mbiot * sigma0 / (Moed2 + p2.alpha^2 * p2.Mbiot);

fprintf('\n=== Table 2 ===\n')
fprintf('  Moed        = %.4g Pa\n', Moed2)
fprintf('  alpha^2*M/Kdr = %.4f  (>>1 => diverges)\n', ratio2)
fprintf('  cv          = %.4g m^2/s\n', cv2)
fprintf('  p0 (undrained) = %.4g Pa\n', p0_2)

%% ---- Build global matrices (Table 2) ----
[K2, H2, Q2, S2] = buildGlobalMatrices(mesh, eleType, p2, npg);
F2 = assembleLoad(mesh, eleType, sigma0, npg);
[freeDofU2, freeDofP2] = buildCaseBC(mesh, W, L);
[U0_2, P0_2] = undrainedIC(K2, Q2, F2, freeDofU2, freeDofP2, p0_2, nNod);

%% ---- Task 2: Staggered Table 2, dt=0.1s, tmax=600s ----
fprintf('\n--- Task 2: staggered Table 2, dt=0.1s ---\n')
dt2    = 0.1;
tmax2  = 600.0;
saveTimes2 = [1, 5, 10, 30, 60, 300, 600];

[U2, P2, t2, iter2, pNorm2] = staggeredSolver(...
    K2, H2, Q2, S2, F2, U0_2, P0_2, freeDofU2, freeDofP2, ...
    dt2, tmax2, tol, maxIter, saveTimes2);

% Plot pressure norm vs iteration for first diverging step
plotPressureNorm(pNorm2, dt2, 'Table 2 - Staggered (divergence)')

%% ---- Task 3: Fixed-stress split, Table 2 ----
fprintf('\n--- Task 3: fixed-stress split, Table 2, dt=0.1s ---\n')
% Contraction factor rho = alpha^2*M/(Kdr + alpha^2*M) ~ 0.89 => ~160 iters for 1e-8
% Use tol=1e-6 => ~120 iterations suffice
Smass2  = p2.Mbiot * S2;                        % int N'N dOmega
Lstab2  = (p2.alpha^2 / Kdr2) * Smass2;         % stabilization matrix
tol3    = 1e-6;
maxIter3 = 200;

[U3, P3, t3, iter3] = fixedStressSolver(...
    K2, H2, Q2, S2, F2, U0_2, P0_2, freeDofU2, freeDofP2, ...
    dt2, tmax2, tol3, maxIter3, saveTimes2, Lstab2);

fprintf('  Mean iterations/step: %.2f\n', mean(iter3(iter3>0)))
plotConsolidation(mesh, U3, P3, t3, p2, 'Table 2 - Fixed-Stress Split')

fprintf('\nDone.\n')
