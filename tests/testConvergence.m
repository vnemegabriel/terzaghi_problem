function tests = testConvergence
% Numerical-accuracy verification: mesh-h and time-step convergence.
tests = functiontests(localfunctions);
end

function p = table1Params()
p.lambda = 4.0e7;
p.mu     = 4.0e7;
p.alpha  = 0.4;
p.Mbiot  = 4.1667e7;
p.kperm  = 1.01937e-9;
p.muf    = 1.0;
p.sigma0 = 10000;
p.L      = 10.0;
end

function err = solveAndError(W, L, nx, ny, eleType, dt, tmax, params)
mesh = generateColumnMesh(W, L, nx, ny, eleType);
[K, H, Q, S] = buildGlobalMatrices(mesh, [], params, 2);
F            = assembleLoad(mesh, eleType, params.sigma0, 2);
[freeU, freeP] = buildCaseBC(mesh, W, L);

Moed = params.lambda + 2*params.mu;
p0   = params.alpha * params.Mbiot * params.sigma0 / (Moed + params.alpha^2 * params.Mbiot);
[U0, P0] = undrainedIC(K, Q, F, freeU, freeP, p0);

[~, P_hist, t_hist, ~, ~] = staggeredSolver( ...
    K, H, Q, S, F, U0, P0, freeU, freeP, ...
    dt, tmax, 1e-10, 100, tmax);

% L2 error along centreline. Pressure DOFs live on corner nodes only,
% so the centreline must be selected within the pressure (corner) space.
cornerNodes = pressureNodeMap(mesh);
xC = mesh.nodes(cornerNodes, 1);
yC = mesh.nodes(cornerNodes, 2);
on_centre = abs(xC) < 1e-10;
yC   = yC(on_centre);
pFEM = P_hist(on_centre, 1);
[ySorted, iS] = sort(yC);
pFEM = pFEM(iS);

x_q = linspace(0, L, 200)';
[p_anal, ~] = analyticalTerzaghi(x_q, t_hist(1), params);
err = sqrt(trapz(ySorted, (pFEM - interp1(x_q, p_anal, ySorted)).^2)) / p0;
end

% ---- mesh-h convergence -----------------------------------------------

function testMeshConvergenceQ4(testCase)
p = table1Params;
nyList = [10 20 40];
err = zeros(numel(nyList), 1);
for i = 1:numel(nyList)
    err(i) = solveAndError(1, p.L, 1, nyList(i), 'Q4', 1, 300, p);
end
% Q4 should be O(h^2) -> slope ~ -2 on log-log of (err vs 1/ny)
h = 1 ./ nyList(:);
loglog_slope = polyfit(log(h), log(err), 1);
slope = loglog_slope(1);
fprintf('Q4 mesh-h slope: %.2f (expected ~2)\n', slope)
verifyGreaterThan(testCase, slope, 1.5)
verifyLessThan(testCase, slope, 2.6)
end

function testMeshConvergenceQ8(testCase)
% Use very small dt so spatial error dominates over time error.
p = table1Params;
nyList = [4 8 16];
err = zeros(numel(nyList), 1);
for i = 1:numel(nyList)
    err(i) = solveAndError(1, p.L, 1, nyList(i), 'Q8', 0.1, 30, p);
end
h = 1 ./ nyList(:);
loglog_slope = polyfit(log(h), log(err), 1);
slope = loglog_slope(1);
fprintf('Q8 mesh-h slope: %.2f (expected >= 2)\n', slope)
% Q8 displacement with Q4 (corner) pressure — Taylor-Hood mixed element.
% Pressure is bilinear, so the pressure L2 rate is ~order 2 here.
verifyGreaterThan(testCase, slope, 1.5)
end

% ---- time-step convergence --------------------------------------------

function testTimeStepConvergenceQ4(testCase)
p = table1Params;
dtList = [2.0 1.0 0.5];
err = zeros(numel(dtList), 1);
for i = 1:numel(dtList)
    err(i) = solveAndError(1, p.L, 1, 40, 'Q4', dtList(i), 100, p);
end
% Backward Euler is O(dt)
loglog_slope = polyfit(log(dtList(:)), log(err), 1);
slope = loglog_slope(1);
fprintf('Q4 dt slope: %.2f (expected ~1)\n', slope)
% Spatial error dominates if too coarse; allow wide band
verifyGreaterThan(testCase, slope, 0.3)
end
