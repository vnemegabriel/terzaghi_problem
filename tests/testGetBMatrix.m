function tests = testGetBMatrix
tests = functiontests(localfunctions);
end

function testBMatrixSize(testCase)
% Q4 element 1x1 in physical space, B should be 3x8.
nodes = [0 0; 1 0; 1 1; 0 1];
[B, detJ] = getBMatrix([0 0], nodes, 'Q4');
verifyEqual(testCase, size(B), [3 8])
verifyGreaterThan(testCase, detJ, 0)
end

function testBMatrixUnitSquare(testCase)
% B at centre of unit square: known closed-form values for Q4
% N_i = 1/4 +/- xi/4 +/- eta/4 +/- xi*eta/4
% At (xi,eta)=(0,0): dN/dx = +/- 1/2 (since Jacobian is 1/2*I)... let's verify by FD.
nodes = [0 0; 2 0; 2 2; 0 2];      % 2x2 square, J = I, detJ = 1
[B, detJ] = getBMatrix([0 0], nodes, 'Q4');
verifyEqual(testCase, detJ, 1.0, 'AbsTol', 1e-12)

% At centre of 2x2 square with corners ordered CCW, dN_i/dx = -/+ 1/4
expected_dNdx = [-0.25  0.25  0.25 -0.25];
expected_dNdy = [-0.25 -0.25  0.25  0.25];
verifyEqual(testCase, B(1, 1:2:end), expected_dNdx, 'AbsTol', 1e-12)
verifyEqual(testCase, B(2, 2:2:end), expected_dNdy, 'AbsTol', 1e-12)
end

function testBMatrixRigidBodyTranslation(testCase)
% Rigid translation produces zero strain.
nodes = [0 0; 1.5 0.2; 1.7 1.3; -0.1 1.4];     % arbitrary quad
u_rigid = repmat([0.7; -0.3], 4, 1);            % every node moved (0.7, -0.3)
[B, ~] = getBMatrix([0.3 -0.5], nodes, 'Q4');
strain = B * u_rigid;
verifyEqual(testCase, strain, zeros(3,1), 'AbsTol', 1e-12)
end

function testBMatrixConstantStrain(testCase)
% Linear displacement field u_x = a*x => eps_xx = a, others zero.
nodes = [0 0; 2 0; 2 1.5; 0 1.5];
a = 0.123;
u = zeros(8,1);
for k = 1:4
    u(2*k-1) = a * nodes(k,1);     % u_x = a*x
    u(2*k)   = 0;
end
[B, ~] = getBMatrix([0.4 -0.2], nodes, 'Q4');
strain = B * u;
verifyEqual(testCase, strain(1), a,  'AbsTol', 1e-12)   % eps_xx
verifyEqual(testCase, strain(2), 0,  'AbsTol', 1e-12)   % eps_yy
verifyEqual(testCase, strain(3), 0,  'AbsTol', 1e-12)   % 2*eps_xy
end
