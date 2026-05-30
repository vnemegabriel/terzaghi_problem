function tests = testGenerateColumnMesh
tests = functiontests(localfunctions);
end

function testMeshNodeCount(testCase)
mesh = generateColumnMesh(1.0, 10.0, 1, 20);
verifyEqual(testCase, size(mesh.nodes,1), 2*21)   % (nx+1)*(ny+1)
verifyEqual(testCase, numel(mesh.elements), 20)   % nx*ny
end

function testMeshExtent(testCase)
W = 1.5; H = 7.3;
mesh = generateColumnMesh(W, H, 3, 5);
verifyEqual(testCase, min(mesh.nodes(:,1)), 0,  'AbsTol', 1e-12)
verifyEqual(testCase, max(mesh.nodes(:,1)), W,  'AbsTol', 1e-12)
verifyEqual(testCase, min(mesh.nodes(:,2)), 0,  'AbsTol', 1e-12)
verifyEqual(testCase, max(mesh.nodes(:,2)), H,  'AbsTol', 1e-12)
end

function testMeshElementCCW(testCase)
% Every Q4 element must have positive area (CCW orientation).
mesh = generateColumnMesh(1.0, 10.0, 2, 8);
for k = 1:numel(mesh.elements)
    n = mesh.nodes(mesh.elements{k}, :);
    % Area = 0.5 * |x1(y2-y4) + x2(y3-y1) + x3(y4-y2) + x4(y1-y3)|
    A = 0.5 * ( n(1,1)*(n(2,2)-n(4,2)) + n(2,1)*(n(3,2)-n(1,2)) ...
              + n(3,1)*(n(4,2)-n(2,2)) + n(4,1)*(n(1,2)-n(3,2)) );
    verifyGreaterThan(testCase, A, 0)
end
end

function testMeshConnectivityValid(testCase)
% All connectivity indices must lie within the node count.
mesh = generateColumnMesh(1.0, 10.0, 1, 20);
nNod = size(mesh.nodes,1);
for k = 1:numel(mesh.elements)
    e = mesh.elements{k};
    verifyEqual(testCase, numel(e), 4)
    verifyTrue(testCase, all(e >= 1 & e <= nNod))
    verifyEqual(testCase, numel(unique(e)), 4)        % no repeated nodes
end
end
