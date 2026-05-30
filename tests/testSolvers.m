function tests = testSolvers
tests = functiontests(localfunctions);
end

% ---- shared problem setup ---------------------------------------------

function s = setupTable1Problem()
% Build a small Table 1 problem (Q4, 1 x 10 elements).
W = 1.0; L = 10.0;
mesh = generateColumnMesh(W, L, 1, 10);
nNod = size(mesh.nodes,1);
nElem = numel(mesh.elements);

s.params.lambda = 4.0e7;
s.params.mu     = 4.0e7;
s.params.alpha  = 0.4;
s.params.Mbiot  = 4.1667e7;
s.params.kperm  = 1.01937e-9;
s.params.muf    = 1.0;
s.params.sigma0 = 10000;
s.params.L      = L;

lam = s.params.lambda; mu = s.params.mu;
C = [lam+2*mu lam 0; lam lam+2*mu 0; 0 0 mu];

nDofU = 2*nNod; nDofP = nNod;
K = sparse(nDofU,nDofU); H = sparse(nDofP,nDofP);
Q = sparse(nDofU,nDofP); S = sparse(nDofP,nDofP);
for k = 1:nElem
    gn = mesh.elements{k};
    ln = mesh.nodes(gn,:);
    nne = numel(gn);
    Ke = assembleMechanical(ln, 'Q4', C, 2);
    [He, Qe, Se] = assemblePoroelastic(ln, 'Q4', ...
        s.params.kperm, s.params.muf, s.params.alpha, s.params.Mbiot, 2);
    dofU = zeros(1, 2*nne);
    dofU(1:2:end) = 2*gn - 1;
    dofU(2:2:end) = 2*gn;
    dofP = gn;
    K(dofU,dofU) = K(dofU,dofU) + Ke;
    H(dofP,dofP) = H(dofP,dofP) + He;
    Q(dofU,dofP) = Q(dofU,dofP) + Qe;
    S(dofP,dofP) = S(dofP,dofP) + Se;
end
F = assembleLoad(mesh, 'Q4', s.params.sigma0, 2);

% Boundary conditions
tol = 1e-10;
fixedU = false(nNod, 2);
for n = 1:nNod
    x = mesh.nodes(n,1); y = mesh.nodes(n,2);
    if abs(x) < tol || abs(x - W) < tol,  fixedU(n,1) = true; end
    if abs(y) < tol,                       fixedU(n,2) = true; end
end
fixedP = false(nNod,1);
for n = 1:nNod
    if abs(mesh.nodes(n,2) - L) < tol,  fixedP(n) = true; end
end
allDofU = 1:2*nNod;  allDofP = 1:nNod;
freeDofU = allDofU(~reshape(fixedU',[],1));
freeDofP = allDofP(~fixedP);

% Initial conditions (undrained)
Moed = lam + 2*mu;
p0   = s.params.alpha * s.params.Mbiot * s.params.sigma0 / (Moed + s.params.alpha^2 * s.params.Mbiot);
P0 = zeros(nNod,1); P0(freeDofP) = p0;
U0 = zeros(2*nNod,1);
rhs = F + Q*P0;
U0(freeDofU) = K(freeDofU,freeDofU) \ rhs(freeDofU);

s.mesh = mesh;
s.K = K; s.H = H; s.Q = Q; s.S = S; s.F = F;
s.U0 = U0; s.P0 = P0;
s.freeDofU = freeDofU; s.freeDofP = freeDofP;
s.p0 = p0;
end

% ---- staggered solver --------------------------------------------------

function testStaggeredConvergesTable1(testCase)
s = setupTable1Problem;
[U_hist, P_hist, t_hist, iter, ~] = staggeredSolver( ...
    s.K, s.H, s.Q, s.S, s.F, s.U0, s.P0, s.freeDofU, s.freeDofP, ...
    1.0, 100.0, 1e-8, 50, [50, 100]);

verifyEqual(testCase, t_hist, [50, 100])
verifyTrue(testCase, all(iter(iter>0) < 50), 'every step should converge')
verifyTrue(testCase, all(isfinite(P_hist(:))))
verifyTrue(testCase, all(isfinite(U_hist(:))))
end

function testStaggeredMatchesAnalyticalTable1(testCase)
s = setupTable1Problem;
[~, P_hist, t_hist, ~, ~] = staggeredSolver( ...
    s.K, s.H, s.Q, s.S, s.F, s.U0, s.P0, s.freeDofU, s.freeDofP, ...
    1.0, 300.0, 1e-8, 50, 300);

yN = s.mesh.nodes(:,2); [yS, iS] = sort(yN);
pFEM = P_hist(iS, 1);
x_q = linspace(0, s.params.L, 200)';
[p_anal, ~] = analyticalTerzaghi(x_q, t_hist(1), s.params);
err = max(abs(pFEM - interp1(x_q, p_anal, yS)));
verifyLessThan(testCase, err / s.p0, 0.01)   % <1% rel error
end

function testStaggeredDivergesTable2(testCase)
% Build same mesh + Table 2; expect divergence within a few steps.
s = setupTable1Problem;     % reuse mesh & DOFs
% Override material to Table 2
s.params.Mbiot = 6.06e9;

% Rebuild only the matrices that depend on Mbiot (S) - regenerate all for clarity
mesh = s.mesh; nNod = size(mesh.nodes,1);
nDofP = nNod;
H = sparse(nDofP,nDofP); Q = sparse(2*nNod,nDofP); S = sparse(nDofP,nDofP);
for k = 1:numel(mesh.elements)
    gn = mesh.elements{k};
    ln = mesh.nodes(gn,:);
    nne = numel(gn);
    [He, Qe, Se] = assemblePoroelastic(ln, 'Q4', ...
        s.params.kperm, s.params.muf, s.params.alpha, s.params.Mbiot, 2);
    dofU = zeros(1, 2*nne);
    dofU(1:2:end) = 2*gn - 1;  dofU(2:2:end) = 2*gn;
    dofP = gn;
    H(dofP,dofP) = H(dofP,dofP) + He;
    Q(dofU,dofP) = Q(dofU,dofP) + Qe;
    S(dofP,dofP) = S(dofP,dofP) + Se;
end

Moed = s.params.lambda + 2*s.params.mu;
p0 = s.params.alpha * s.params.Mbiot * s.params.sigma0 / (Moed + s.params.alpha^2 * s.params.Mbiot);
P0 = zeros(nNod,1); P0(s.freeDofP) = p0;
U0 = zeros(2*nNod,1);
rhs = s.F + Q*P0;
U0(s.freeDofU) = s.K(s.freeDofU,s.freeDofU) \ rhs(s.freeDofU);

[~, P_hist, t_hist, ~, pNorm] = staggeredSolver( ...
    s.K, H, Q, S, s.F, U0, P0, s.freeDofU, s.freeDofP, ...
    0.1, 5.0, 1e-8, 50, 0.5:0.5:5.0);

% Either: solution blew up and t_hist is truncated,
% OR: pressure norms grew across iterations in step 1.
hasDiverged = numel(t_hist) < 10;
nrm1 = pNorm{1};
growing = numel(nrm1) > 5 && nrm1(end) > 10 * nrm1(1);
verifyTrue(testCase, hasDiverged || growing, ...
    'Table 2 staggered should either blow up or show growing pressure norms')
end

% ---- fixed-stress solver -----------------------------------------------

function testFixedStressMatchesMonolithic(testCase)
s = setupTable1Problem;        % use Table 1 (easy case)
% Build Lstab
Kdr = s.params.lambda + 2*s.params.mu;
Smass = s.params.Mbiot * s.S;
Lstab = (s.params.alpha^2 / Kdr) * Smass;

dt = 1.0; tmax = 60;
[~, P_fs, t_fs, ~] = fixedStressSolver( ...
    s.K, s.H, s.Q, s.S, s.F, s.U0, s.P0, s.freeDofU, s.freeDofP, ...
    dt, tmax, 1e-10, 100, 60, Lstab);

% Monolithic reference
nDofU = numel(s.U0); nDofP = numel(s.P0);
nTot = nDofU + nDofP;
A_mono = [s.K, -s.Q; s.Q'/dt, s.S/dt + s.H];
freeAllU = false(nDofU,1); freeAllU(s.freeDofU) = true;
freeAllP = false(nDofP,1); freeAllP(s.freeDofP) = true;
freeAll = [freeAllU; freeAllP];

Un = s.U0; Pn = s.P0;
for step = 1:60
    rhs = [s.F; s.S/dt*Pn + s.Q'/dt*Un];
    sol = zeros(nTot,1);
    sol(freeAll) = A_mono(freeAll,freeAll) \ rhs(freeAll);
    Un = sol(1:nDofU);
    Pn = sol(nDofU+1:end);
end
err = max(abs(P_fs(:,1) - Pn));
verifyLessThan(testCase, err / s.p0, 1e-6)
end

function testFixedStressMatchesAnalyticalTable2(testCase)
% Set up Table 2 problem and check fixed-stress matches analytical.
s = setupTable1Problem;
s.params.Mbiot = 6.06e9;

mesh = s.mesh; nNod = size(mesh.nodes,1);
H = sparse(nNod,nNod); Q = sparse(2*nNod,nNod); S = sparse(nNod,nNod);
for k = 1:numel(mesh.elements)
    gn = mesh.elements{k};
    ln = mesh.nodes(gn,:);
    nne = numel(gn);
    [He, Qe, Se] = assemblePoroelastic(ln, 'Q4', ...
        s.params.kperm, s.params.muf, s.params.alpha, s.params.Mbiot, 2);
    dofU = zeros(1, 2*nne);
    dofU(1:2:end) = 2*gn - 1; dofU(2:2:end) = 2*gn;
    dofP = gn;
    H(dofP,dofP) = H(dofP,dofP) + He;
    Q(dofU,dofP) = Q(dofU,dofP) + Qe;
    S(dofP,dofP) = S(dofP,dofP) + Se;
end

Moed = s.params.lambda + 2*s.params.mu;
Kdr  = Moed;
p0 = s.params.alpha * s.params.Mbiot * s.params.sigma0 / (Moed + s.params.alpha^2 * s.params.Mbiot);
P0 = zeros(nNod,1); P0(s.freeDofP) = p0;
U0 = zeros(2*nNod,1);
rhs = s.F + Q*P0;
U0(s.freeDofU) = s.K(s.freeDofU,s.freeDofU) \ rhs(s.freeDofU);

Lstab = (s.params.alpha^2 / Kdr) * (s.params.Mbiot * S);

[~, P_hist, t_hist, ~] = fixedStressSolver( ...
    s.K, H, Q, S, s.F, U0, P0, s.freeDofU, s.freeDofP, ...
    0.1, 300, 1e-8, 200, 300, Lstab);

yN = mesh.nodes(:,2); [yS, iS] = sort(yN);
pFEM = P_hist(iS, 1);
x_q = linspace(0, s.params.L, 200)';
[p_anal, ~] = analyticalTerzaghi(x_q, t_hist(1), s.params);
err = max(abs(pFEM - interp1(x_q, p_anal, yS)));
verifyLessThan(testCase, err / p0, 0.02)   % <2% rel error
end
