function tests = testFuzz
% Randomised stress tests: generate N configurations and check invariants.
% Uses a fixed RNG seed for reproducibility.
tests = functiontests(localfunctions);
end

function testFuzzAssemblyInvariants(testCase)
% For each random configuration check:
%   * K, H, S are symmetric (to fp tolerance)
%   * eig(K(free,free)) > 0
%   * sum(F_y) == -sigma0 * W

rng(42)
nRuns = 20;

failures = {};
for r = 1:nRuns
    cfg = randomConfig();
    try
        mesh = generateColumnMesh(cfg.W, cfg.L, cfg.nx, cfg.ny, cfg.eleType);
        [K, H, ~, S] = buildGlobalMatrices(mesh, [], cfg.params, 2);
        F = assembleLoad(mesh, cfg.eleType, cfg.sigma0, 2);
        [freeU, ~] = buildCaseBC(mesh, cfg.W, cfg.L);

        % Symmetry
        if norm(K - K', 'fro') > 1e-4 * norm(K,'fro')
            failures{end+1} = sprintf('run %d: K not symmetric', r); %#ok<AGROW>
        end
        if norm(H - H', 'fro') > 1e-10
            failures{end+1} = sprintf('run %d: H not symmetric', r); %#ok<AGROW>
        end
        if norm(S - S', 'fro') > 1e-10
            failures{end+1} = sprintf('run %d: S not symmetric', r); %#ok<AGROW>
        end

        % Positive-definiteness of K on free DOFs
        Kf = K(freeU, freeU);
        [~, pp] = chol(Kf);
        if pp ~= 0
            failures{end+1} = sprintf('run %d: K(free) not SPD', r); %#ok<AGROW>
        end

        % Load conservation
        F_y_total = sum(F(2:2:end));
        if abs(F_y_total + cfg.sigma0 * cfg.W) > 1e-8 * cfg.sigma0 * cfg.W
            failures{end+1} = sprintf('run %d: load not conserved', r); %#ok<AGROW>
        end
    catch ME
        failures{end+1} = sprintf('run %d: unexpected error %s', r, ME.identifier); %#ok<AGROW>
    end
end

if ~isempty(failures)
    fprintf('Fuzz assembly failures:\n')
    for i = 1:numel(failures)
        fprintf('  %s\n', failures{i})
    end
end
verifyEmpty(testCase, failures)
end

function testFuzzAnalyticalBounded(testCase)
% Analytical pressure must be bounded for arbitrary inputs.
rng(42)
nRuns = 20;

failures = {};
for r = 1:nRuns
    cfg = randomConfig();
    p   = cfg.params; p.sigma0 = cfg.sigma0; p.L = cfg.L;

    x = linspace(0, cfg.L, 30)';
    t = 10.^(rand*3);                        % 1 .. 1000 s
    try
        [pa, ~] = analyticalTerzaghi(x, t, p);
        Moed = p.lambda + 2*p.mu;
        p0   = p.alpha*p.Mbiot*p.sigma0 / (Moed + p.alpha^2*p.Mbiot);
        if any(pa < -0.2*p0) || any(pa > 1.2*p0)
            failures{end+1} = sprintf('run %d: bound violation (range [%.2g, %.2g], p0=%.2g)', ...
                r, min(pa), max(pa), p0); %#ok<AGROW>
        end
        if any(~isfinite(pa))
            failures{end+1} = sprintf('run %d: NaN/Inf in pressure', r); %#ok<AGROW>
        end
    catch ME
        failures{end+1} = sprintf('run %d: unexpected error %s', r, ME.identifier); %#ok<AGROW>
    end
end

if ~isempty(failures)
    fprintf('Fuzz analytical failures:\n')
    for i = 1:numel(failures), fprintf('  %s\n', failures{i}); end
end
verifyEmpty(testCase, failures)
end

% ---- helpers -----------------------------------------------------------

function cfg = randomConfig()
cfg.nx = randi([1 4]);
cfg.ny = randi([5 15]);
cfg.W  = 1 + 9*rand;
cfg.L  = 1 + 9*rand;
if rand > 0.5, cfg.eleType = 'Q4'; else, cfg.eleType = 'Q8'; end

cfg.params.lambda = 10^(6 + 2*rand);
cfg.params.mu     = 10^(6 + 2*rand);
cfg.params.alpha  = 0.1 + 0.8*rand;
cfg.params.Mbiot  = 10^(7 + 2*rand);
cfg.params.kperm  = 10^(-12 + 4*rand);
cfg.params.muf    = 0.5 + rand;
cfg.sigma0        = 1000 + 1e4*rand;
end
