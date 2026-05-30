function tests = testQ8Mesh
% Q8-specific tests: mesh generation, element-count, connectivity, and
% end-to-end assembly + solve.
tests = functiontests(localfunctions);
end

% ---- mesh structure ----------------------------------------------------

function testQ8NodeCount(testCase)
% nNod = (ny+1)*(2*nx+1) + ny*(nx+1)
mesh = generateColumnMesh(1.0, 10.0, 1, 20, 'Q8');
expected = 21*3 + 20*2;          % 63 + 40
verifyEqual(testCase, size(mesh.nodes,1), expected)
verifyEqual(testCase, numel(mesh.elements), 20)
end

function testQ8ConnectivitySize(testCase)
mesh = generateColumnMesh(1.0, 10.0, 2, 5, 'Q8');
for k = 1:numel(mesh.elements)
    verifyEqual(testCase, numel(mesh.elements{k}), 8)
    verifyEqual(testCase, numel(unique(mesh.elements{k})), 8)
end
end

function testQ8MidEdgeNodesOnEdges(testCase)
% Each mid-edge node must lie exactly halfway between its two corner nodes.
mesh = generateColumnMesh(1.0, 10.0, 2, 3, 'Q8');
for k = 1:numel(mesh.elements)
    e = mesh.elements{k};
    p = mesh.nodes(e, :);                  % [8 x 2]
    % n5 between n1 and n2 (bottom)
    verifyEqual(testCase, p(5,:), 0.5*(p(1,:)+p(2,:)), 'AbsTol', 1e-12)
    % n6 between n2 and n3 (right)
    verifyEqual(testCase, p(6,:), 0.5*(p(2,:)+p(3,:)), 'AbsTol', 1e-12)
    % n7 between n3 and n4 (top)
    verifyEqual(testCase, p(7,:), 0.5*(p(3,:)+p(4,:)), 'AbsTol', 1e-12)
    % n8 between n4 and n1 (left)
    verifyEqual(testCase, p(8,:), 0.5*(p(4,:)+p(1,:)), 'AbsTol', 1e-12)
end
end

function testQ8SharedMidEdgesBetweenElements(testCase)
% Adjacent vertical elements must share the horizontal mid-edge nodes
% (top of lower element == bottom of upper element).
mesh = generateColumnMesh(1.0, 10.0, 1, 3, 'Q8');
for j = 1:2
    lower = mesh.elements{j};
    upper = mesh.elements{j+1};
    % Corners shared: lower n3 == upper n2, lower n4 == upper n1
    verifyEqual(testCase, lower(3), upper(2))
    verifyEqual(testCase, lower(4), upper(1))
    % Horizontal mid-edge: lower n7 (top mid) == upper n5 (bottom mid)
    verifyEqual(testCase, lower(7), upper(5))
end
end

% ---- assembly with Q8 --------------------------------------------------

function testQ8KeSizeAndSymmetry(testCase)
mesh = generateColumnMesh(1.0, 1.0, 1, 1, 'Q8');
nodes = mesh.nodes(mesh.elements{1}, :);
lam = 4e7; mu = 4e7;
C = [lam+2*mu lam 0; lam lam+2*mu 0; 0 0 mu];
Ke = assembleMechanical(nodes, 'Q8', C, 2);
verifyEqual(testCase, size(Ke), [16 16])
verifyLessThan(testCase, norm(Ke - Ke', 'fro') / norm(Ke,'fro'), 1e-10)
end

function testQ8PoroelasticSizes(testCase)
mesh = generateColumnMesh(1.0, 1.0, 1, 1, 'Q8');
nodes = mesh.nodes(mesh.elements{1}, :);
[He, Qe, Se] = assemblePoroelastic(nodes, 'Q8', 1e-9, 1.0, 0.4, 4.17e7, 2);
verifyEqual(testCase, size(He), [8  8])
verifyEqual(testCase, size(Qe), [16 8])
verifyEqual(testCase, size(Se), [8  8])
verifyLessThan(testCase, norm(He - He', 'fro'), 1e-15)
verifyLessThan(testCase, norm(Se - Se', 'fro'), 1e-15)
end

function testQ8LoadTotalForce(testCase)
% Total load on top face must equal -sigma0 * W regardless of element type.
mesh = generateColumnMesh(1.0, 10.0, 2, 5, 'Q8');
sigma0 = 10000;
F = assembleLoad(mesh, 'Q8', sigma0, 2);
F_y_total = sum(F(2:2:end));
verifyEqual(testCase, F_y_total, -sigma0 * 1.0, 'RelTol', 1e-12)
end

% ---- end-to-end solve --------------------------------------------------

function testQ8FullPipelineMatchesAnalytical(testCase)
% Full Q8 staggered solve must match analytical solution within 1% on Table 1.
W = 1.0; L = 10.0;
mesh = generateColumnMesh(W, L, 1, 10, 'Q8');

p.lambda = 4.0e7;
p.mu     = 4.0e7;
p.alpha  = 0.4;
p.Mbiot  = 4.1667e7;
p.kperm  = 1.01937e-9;
p.muf    = 1.0;
p.sigma0 = 10000;
p.L      = L;

[K, H, Q, S] = buildGlobalMatrices(mesh, 'Q8', p, 2);
F            = assembleLoad(mesh, 'Q8', p.sigma0, 2);
[freeU, freeP] = buildCaseBC(mesh, W, L);

Moed = p.lambda + 2*p.mu;
p0   = p.alpha * p.Mbiot * p.sigma0 / (Moed + p.alpha^2 * p.Mbiot);
nNod = size(mesh.nodes,1);
[U0, P0] = undrainedIC(K, Q, F, freeU, freeP, p0, nNod);

[~, P_hist, t_hist, iter, ~] = staggeredSolver( ...
    K, H, Q, S, F, U0, P0, freeU, freeP, ...
    1.0, 300.0, 1e-8, 50, 300);

verifyTrue(testCase, all(iter(iter>0) < 50))

% Compare centreline pressures
on_centre = abs(mesh.nodes(:,1)) < 1e-10;
yC = mesh.nodes(on_centre, 2);
pFEM = P_hist(on_centre, 1);
[ySorted, iS] = sort(yC);
pFEM = pFEM(iS);

x_q = linspace(0, L, 200)';
[p_anal, ~] = analyticalTerzaghi(x_q, t_hist(1), p);
err = max(abs(pFEM - interp1(x_q, p_anal, ySorted)));
verifyLessThan(testCase, err / p0, 0.01)
end
