function tests = testShapeFunctions
% Test suite for shape function library.
tests = functiontests(localfunctions);
end

% ---- Q4 ----------------------------------------------------------------

function testQ4PartitionOfUnity(testCase)
% Sum of all shape functions = 1 anywhere.
pts = [-0.7 0.3; 0 0; 0.5 -0.8; 1 1; -1 -1];
for i = 1:size(pts,1)
    N = shapeFunctions(pts(i,:), 'Q4');
    verifyEqual(testCase, sum(N), 1.0, 'AbsTol', 1e-12)
end
end

function testQ4DeltaProperty(testCase)
% N_i is 1 at node i, 0 at other corners.
corners = [-1 -1; 1 -1; 1 1; -1 1];
for i = 1:4
    N = shapeFunctions(corners(i,:), 'Q4');
    expected = zeros(1,4);  expected(i) = 1;
    verifyEqual(testCase, N, expected, 'AbsTol', 1e-12)
end
end

function testQ4DerivativeSumZero(testCase)
% Sum of derivatives must be zero (partition of unity is constant).
pts = [-0.5 0.5; 0.7 -0.2; 0 0];
for i = 1:size(pts,1)
    dN = shapeFunctionsDer(pts(i,:), 'Q4');
    verifyEqual(testCase, sum(dN(1,:)), 0, 'AbsTol', 1e-12)
    verifyEqual(testCase, sum(dN(2,:)), 0, 'AbsTol', 1e-12)
end
end

function testQ4DerivativesByFiniteDifference(testCase)
% Spot-check dN/dxi by central difference.
h = 1e-6;
xi = 0.3; eta = -0.4;
dN_analytic = shapeFunctionsDer([xi eta], 'Q4');
Np = shapeFunctions([xi+h eta], 'Q4');
Nm = shapeFunctions([xi-h eta], 'Q4');
dN_dxi_fd = (Np - Nm) / (2*h);
verifyEqual(testCase, dN_analytic(1,:), dN_dxi_fd, 'AbsTol', 1e-8)

Np = shapeFunctions([xi eta+h], 'Q4');
Nm = shapeFunctions([xi eta-h], 'Q4');
dN_deta_fd = (Np - Nm) / (2*h);
verifyEqual(testCase, dN_analytic(2,:), dN_deta_fd, 'AbsTol', 1e-8)
end

% ---- Q8 ----------------------------------------------------------------

function testQ8PartitionOfUnity(testCase)
pts = [-0.7 0.3; 0 0; 0.5 -0.8];
for i = 1:size(pts,1)
    N = shapeFunctions(pts(i,:), 'Q8');
    verifyEqual(testCase, sum(N), 1.0, 'AbsTol', 1e-12)
end
end

function testQ8DeltaProperty(testCase)
% Corner + mid-edge nodes for Q8.
nodes_iso = [-1 -1;  1 -1;  1  1; -1  1; ...
              0 -1;  1  0;  0  1; -1  0];
for i = 1:8
    N = shapeFunctions(nodes_iso(i,:), 'Q8');
    expected = zeros(1,8);  expected(i) = 1;
    verifyEqual(testCase, N, expected, 'AbsTol', 1e-12)
end
end

% ---- CST ---------------------------------------------------------------

function testCSTPartitionOfUnity(testCase)
% Inside reference triangle: xi+eta <= 1, both >= 0
pts = [0.1 0.2; 0.3 0.5; 0.0 0.0; 1.0 0.0; 0.0 1.0];
for i = 1:size(pts,1)
    N = shapeFunctions(pts(i,:), 'CST');
    verifyEqual(testCase, sum(N), 1.0, 'AbsTol', 1e-12)
end
end
