function tests = testAssembly
tests = functiontests(localfunctions);
end

function testKeSymmetric(testCase)
nodes = [0 0; 1 0; 1 1; 0 1];
lam = 4e7; mu = 4e7;
C = [lam+2*mu lam 0; lam lam+2*mu 0; 0 0 mu];
Ke = assembleMechanical(nodes, 'Q4', C, 2);
verifyLessThan(testCase, norm(Ke - Ke', 'fro'), 1e-6)
end

function testKeRigidBodyNullSpace(testCase)
% Rigid-body translations must produce zero internal force.
nodes = [0 0; 2 0.3; 2.1 1.7; -0.1 1.8];
lam = 4e7; mu = 4e7;
C = [lam+2*mu lam 0; lam lam+2*mu 0; 0 0 mu];
Ke = assembleMechanical(nodes, 'Q4', C, 2);

u_tx = repmat([1; 0], 4, 1);
u_ty = repmat([0; 1], 4, 1);
verifyLessThan(testCase, norm(Ke * u_tx), 1e-6)
verifyLessThan(testCase, norm(Ke * u_ty), 1e-6)
end

function testKeConstantStrainEnergy(testCase)
% For a unit square loaded by pure xx strain a, energy = 0.5 * (lam+2mu) * a^2 * Area
nodes = [0 0; 1 0; 1 1; 0 1];
lam = 4e7; mu = 4e7;
C = [lam+2*mu lam 0; lam lam+2*mu 0; 0 0 mu];
Ke = assembleMechanical(nodes, 'Q4', C, 2);

a = 1e-4;
u = zeros(8,1);
for k = 1:4
    u(2*k-1) = a * nodes(k,1);
end
energy_num   = 0.5 * u' * Ke * u;
energy_exact = 0.5 * (lam + 2*mu) * a^2 * 1.0;     % area = 1
verifyEqual(testCase, energy_num, energy_exact, 'RelTol', 1e-10)
end

function testHeSymmetric(testCase)
nodes = [0 0; 1 0; 1 1; 0 1];
[He, ~, ~] = assemblePoroelastic(nodes, 'Q4', 1e-9, 1.0, 0.4, 4.17e7, 2);
verifyLessThan(testCase, norm(He - He', 'fro'), 1e-15)
end

function testSeSymmetricPositive(testCase)
nodes = [0 0; 1 0; 1 1; 0 1];
[~, ~, Se] = assemblePoroelastic(nodes, 'Q4', 1e-9, 1.0, 0.4, 4.17e7, 2);
verifyLessThan(testCase, norm(Se - Se', 'fro'), 1e-15)
verifyGreaterThan(testCase, min(eig(Se)), 0)
end

function testQeSize(testCase)
nodes = [0 0; 1 0; 1 1; 0 1];
[~, Qe, ~] = assemblePoroelastic(nodes, 'Q4', 1e-9, 1.0, 0.4, 4.17e7, 2);
verifyEqual(testCase, size(Qe), [8 4])
end

function testLoadAssemblyTotalForce(testCase)
% Total y-force on top face must equal -sigma0 * W
W = 1.0; L = 10.0;
sigma0 = 10000;
mesh = generateColumnMesh(W, L, 2, 5);
F = assembleLoad(mesh, 'Q4', sigma0, 2);
F_y_total = sum(F(2:2:end));
verifyEqual(testCase, F_y_total, -sigma0 * W, 'RelTol', 1e-12)
% No x-force expected
F_x_total = sum(F(1:2:end));
verifyEqual(testCase, F_x_total, 0, 'AbsTol', 1e-12)
end

function testLoadOnlyTopFace(testCase)
% Only nodes on the top face should receive a non-zero load.
W = 1.0; L = 10.0;
mesh = generateColumnMesh(W, L, 2, 5);
F = assembleLoad(mesh, 'Q4', 1000, 2);
F_mat = reshape(F, 2, [])';      % [nNod x 2]
tol_node = 1e-10;
for n = 1:size(mesh.nodes,1)
    if mesh.nodes(n,2) < L - tol_node
        verifyEqual(testCase, F_mat(n,:), [0 0], 'AbsTol', 1e-12)
    end
end
end
