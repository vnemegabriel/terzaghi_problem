function tests = testAnalyticalTerzaghi
tests = functiontests(localfunctions);
end

function p = defaultParams()
p.lambda = 4.0e7;
p.mu     = 4.0e7;
p.alpha  = 0.4;
p.Mbiot  = 4.1667e7;
p.kperm  = 1.01937e-9;
p.muf    = 1.0;
p.sigma0 = 10000;
p.L      = 10.0;
end

function testInitialPressureUniform(testCase)
% At t -> 0+, pressure tends to p0 everywhere except at x=L.
% With a truncated 50-term Fourier series, Gibbs oscillations of ~9 %
% appear near x=L. We restrict the check to x <= 0.25*L where the
% series converges well, and use a relaxed tolerance.
p = defaultParams;
x = linspace(0, 0.25*p.L, 20)';
[p_anal, ~] = analyticalTerzaghi(x, 1e-8, p);

Moed = p.lambda + 2*p.mu;
p0   = p.alpha * p.Mbiot * p.sigma0 / (Moed + p.alpha^2 * p.Mbiot);

verifyEqual(testCase, p_anal, p0 * ones(numel(x),1), 'RelTol', 5e-2)
end

function testDrainedBoundaryZero(testCase)
% Pressure at x=L must be 0 for all times.
p = defaultParams;
times = [1, 10, 100, 1000];
p_top = zeros(numel(times),1);
for i = 1:numel(times)
    [pa, ~] = analyticalTerzaghi(p.L, times(i), p);
    p_top(i) = pa;
end
verifyEqual(testCase, p_top, zeros(numel(times),1), 'AbsTol', 1e-9)
end

function testPressureDecaysWithTime(testCase)
% Average pressure must decrease monotonically.
p = defaultParams;
x = linspace(0, p.L, 200)';
times = [10, 100, 1000, 10000];
avg_p = zeros(numel(times),1);
for i = 1:numel(times)
    [pa, ~] = analyticalTerzaghi(x, times(i), p);
    avg_p(i) = mean(pa);
end
verifyTrue(testCase, all(diff(avg_p) < 0))
end

function testPressureBounded(testCase)
% Pressure should always lie between 0 and p0 for the consolidation regime.
p = defaultParams;
Moed = p.lambda + 2*p.mu;
p0   = p.alpha * p.Mbiot * p.sigma0 / (Moed + p.alpha^2 * p.Mbiot);
x = linspace(0, p.L, 50)';
[pa, ~] = analyticalTerzaghi(x, 500, p);
verifyTrue(testCase, all(pa >= -1e-6))
verifyTrue(testCase, all(pa <= p0 + 1e-6))
end

function testFinalPressureZero(testCase)
% At very long time, pressure -> 0 everywhere.
p = defaultParams;
x = linspace(0, p.L, 50)';
[pa, ~] = analyticalTerzaghi(x, 1e7, p);
verifyEqual(testCase, pa, zeros(numel(x),1), 'AbsTol', 1e-3)
end

function testDisplacementAtBottomZero(testCase)
% u(0,t) = 0 by definition (integration from 0).
p = defaultParams;
[~, u_anal] = analyticalTerzaghi(0, 100, p);
verifyEqual(testCase, u_anal, 0, 'AbsTol', 1e-12)
end

function testDisplacementFinalEquilibrium(testCase)
% At t -> infinity, drained settlement: u(x) = -sigma0 * x / Moed.
p = defaultParams;
x = linspace(0, p.L, 20)';
[~, u_inf] = analyticalTerzaghi(x, 1e7, p);
Moed = p.lambda + 2*p.mu;
u_expected = -p.sigma0 * x / Moed;
verifyEqual(testCase, u_inf, u_expected, 'AbsTol', 1e-6)
end
