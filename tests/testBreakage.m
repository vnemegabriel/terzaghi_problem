function tests = testBreakage
% Adversarial / negative tests covering input-validation guards
% added in Phase E3 of the optimization plan.
tests = functiontests(localfunctions);
end

% ---- gauss1D -----------------------------------------------------------

function testGauss1DRejectsZero(testCase)
verifyError(testCase, @() gauss1D(0),  'gauss1D:invalidN')
end

function testGauss1DRejectsTooMany(testCase)
verifyError(testCase, @() gauss1D(7),  'gauss1D:invalidN')
end

% ---- shapeFunctions / shapeFunctionsDer --------------------------------

function testShapeFunctionsRejectsUnknownEleType(testCase)
verifyError(testCase, @() shapeFunctions([0 0], 'XX'), ...
    'shapeFunctions:unknownEleType')
end

function testShapeFunctionsDerRejectsUnknownEleType(testCase)
verifyError(testCase, @() shapeFunctionsDer([0 0], 'P1'), ...
    'shapeFunctionsDer:unknownEleType')
end

function testShapeFunctionsRejectsScalarLocation(testCase)
verifyError(testCase, @() shapeFunctions(0.5, 'Q4'), ...
    'shapeFunctions:badLocation')
end

function testShapeFunctionsDerRejectsScalarLocation(testCase)
verifyError(testCase, @() shapeFunctionsDer(0.5, 'Q4'), ...
    'shapeFunctionsDer:badLocation')
end

% ---- generateColumnMesh ------------------------------------------------

function testMeshRejectsZeroNx(testCase)
verifyError(testCase, @() generateColumnMesh(1, 10, 0, 5), ...
    'MATLAB:generateColumnMesh:expectedPositive')
end

function testMeshRejectsNegativeWidth(testCase)
verifyError(testCase, @() generateColumnMesh(-1, 10, 1, 5), ...
    'MATLAB:generateColumnMesh:expectedPositive')
end

function testMeshRejectsUnknownEleType(testCase)
verifyError(testCase, @() generateColumnMesh(1, 10, 1, 5, 'Q1'), ...
    'generateColumnMesh:invalidEleType')
end

% ---- assembleMechanical ------------------------------------------------

function testAssembleMechanicalRaisesOnTypeMismatch(testCase)
% Q4 mesh nodes (4 rows) processed as Q8 -> shape-fn dim mismatch
nodes = [0 0; 1 0; 1 1; 0 1];
lam = 4e7; mu = 4e7;
C = [lam+2*mu lam 0; lam lam+2*mu 0; 0 0 mu];
verifyError(testCase, ...
    @() assembleMechanical(nodes, 'Q8', C, 2), ...
    'MATLAB:innerdim')
end

% ---- assembleLoad ------------------------------------------------------

function testAssembleLoadRejectsUnknownEleType(testCase)
mesh = generateColumnMesh(1, 10, 1, 5);
verifyError(testCase, @() assembleLoad(mesh, 'XX', 1000, 2), ...
    'assembleLoad:unknownEleType')
end

% ---- solver inputs -----------------------------------------------------

function testStaggeredRejectsZeroDt(testCase)
s = quickProblem();
verifyError(testCase, ...
    @() staggeredSolver(s.K, s.H, s.Q, s.S, s.F, s.U0, s.P0, ...
                        s.freeU, s.freeP, 0, 10, 1e-8, 10, 5), ...
    'MATLAB:staggeredSolver:expectedPositive')
end

function testStaggeredRejectsTmaxBelowDt(testCase)
s = quickProblem();
verifyError(testCase, ...
    @() staggeredSolver(s.K, s.H, s.Q, s.S, s.F, s.U0, s.P0, ...
                        s.freeU, s.freeP, 10, 1, 1e-8, 10, 1), ...
    'staggeredSolver:tmaxBelowDt')
end

function testStaggeredRejectsEmptyFreeDof(testCase)
s = quickProblem();
verifyError(testCase, ...
    @() staggeredSolver(s.K, s.H, s.Q, s.S, s.F, s.U0, s.P0, ...
                        [], s.freeP, 1, 10, 1e-8, 10, 5), ...
    'staggeredSolver:noFreeMechDof')
end

function testFixedStressRejectsBadLstab(testCase)
s = quickProblem();
verifyError(testCase, ...
    @() fixedStressSolver(s.K, s.H, s.Q, s.S, s.F, s.U0, s.P0, ...
                          s.freeU, s.freeP, 1, 10, 1e-8, 10, 5, ...
                          sparse(3,3)), ...
    'fixedStressSolver:LstabSize')
end

% ---- analyticalTerzaghi bounds (positive test) -------------------------

function testAnalyticalIsBounded(testCase)
% Confirms the documented Gibbs ringing is not unbounded.
p.lambda = 4e7; p.mu = 4e7; p.alpha = 0.4; p.Mbiot = 4.17e7;
p.kperm = 1e-9; p.muf = 1; p.sigma0 = 10000; p.L = 10;
x = linspace(0, 10, 50)';
[pa, ~] = analyticalTerzaghi(x, 1e-6, p);
Moed = p.lambda + 2*p.mu;
p0   = p.alpha*p.Mbiot*p.sigma0/(Moed + p.alpha^2*p.Mbiot);
verifyTrue(testCase, all(pa >= -0.2*p0))       % Gibbs lower
verifyTrue(testCase, all(pa <=  1.2*p0))       % Gibbs upper
end

% ---- helpers -----------------------------------------------------------

function s = quickProblem()
% Minimal Q4 problem to feed the solvers.
W = 1; L = 10;
mesh = generateColumnMesh(W, L, 1, 4);
nNod = size(mesh.nodes,1);
lam = 4e7; mu = 4e7;
C = [lam+2*mu lam 0; lam lam+2*mu 0; 0 0 mu];
nDofU = 2*nNod; nDofP = nNod;
K = sparse(nDofU,nDofU); H = sparse(nDofP,nDofP);
Q = sparse(nDofU,nDofP); S = sparse(nDofP,nDofP);
for k = 1:numel(mesh.elements)
    gn = mesh.elements{k};
    ln = mesh.nodes(gn,:);
    nne = numel(gn);
    Ke = assembleMechanical(ln,'Q4',C,2);
    [He,Qe,Se] = assemblePoroelastic(ln,'Q4',1e-9,1,0.4,4.17e7,2);
    dofU = zeros(1,2*nne);
    dofU(1:2:end) = 2*gn-1; dofU(2:2:end) = 2*gn;
    dofP = gn;
    K(dofU,dofU) = K(dofU,dofU)+Ke;
    H(dofP,dofP) = H(dofP,dofP)+He;
    Q(dofU,dofP) = Q(dofU,dofP)+Qe;
    S(dofP,dofP) = S(dofP,dofP)+Se;
end
F = assembleLoad(mesh,'Q4',1000,2);
[freeU, freeP] = buildCaseBC(mesh, W, L);
P0 = zeros(nNod,1); P0(freeP) = 100;
U0 = zeros(2*nNod,1);
rhs = F + Q*P0;
U0(freeU) = K(freeU,freeU)\rhs(freeU);
s.K = K; s.H = H; s.Q = Q; s.S = S; s.F = F;
s.U0 = U0; s.P0 = P0; s.freeU = freeU; s.freeP = freeP;
end
