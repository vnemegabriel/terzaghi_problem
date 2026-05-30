clear;
clc;
close all;

%% 1. MESH DEFINITION

mesh.nodes = 1000 * [
    0 0; 
    0.5 0;
    1 0;
    1.5 0;
    2 0;
    0 0.5; 
    1 0.5;
    2 0.5;
    0 1; 
    0.5 1;
    1 1;
    1.5 1; 
    2 1;
    0 1.5;
    1 1.5; 
    2 1.5;
    0 2; 
    0.5 2; 
    1 2; 
    1.5 2; 
    2 2];

% MISTAKE 1 FIX: Use sequential node ordering to match shapefuns.m
% The node order is now counter-clockwise: corner, midside, corner, etc.
mesh.elements = [
    1,  2,  3,  7, 11, 10,  9,  6;
    3,  4,  5,  8, 13, 12, 11,  7;
    9, 10, 11, 15, 19, 18, 17, 14;
    11, 12, 13, 16, 21, 20, 19, 15];

meshplot(mesh.elements,mesh.nodes,'b',1); % Visualize the mesh

%% 2. GENERAL DEFINITIONS & MATERIAL PROPERTIES
nNodElem = 8; nDofNod = 2;
nNod = size(mesh.nodes, 1);
nDofTot = nDofNod * nNod;
nElem = size(mesh.elements, 1);
nNodElemEdge = 3;

young = 2e5; pnu = 0.33;
t = 1; % Assume unit thickness

E = young / (1 - pnu^2) * [1 pnu 0; pnu 1 0; 0 0 0.5 * (1 - pnu)];
sigma_y_target = 400;

%% 3. Load distro calculations
fprintf('3. Calculating consistent nodal loads for sigma_y traction...\n');
F = zeros(nDofTot, 1);

traction = [0; sigma_y_target]; % Traction in y-direction

top_edge_nodes = [17, 18, 19;
                  19, 20, 21];

[w, gp] = gauss1D(2);

% This loop is correct and calculates the forces accurately.
for i = 1:size(top_edge_nodes, 1)
    edge_nodes_indices = top_edge_nodes(i, :);
    edge_nodes_coords = mesh.nodes(edge_nodes_indices, :);

    Fe = zeros(nNodElemEdge * nDofNod, 1);

    for k = 1:length(gp)
        xi = gp(k);
        % Using standard 1D quadratic shape functions for simplicity and robustness
        N_edge = [0.5*xi*(xi-1), (1-xi^2), 0.5*xi*(xi+1)];
        dN_dxi_edge = [xi-0.5, -2*xi, xi+0.5];

        dx_dxi = dN_dxi_edge * edge_nodes_coords(:, 1);
        dy_dxi = dN_dxi_edge * edge_nodes_coords(:, 2);
        J_edge = sqrt(dx_dxi^2 + dy_dxi^2);

        N_matrix = zeros(2, nNodElemEdge * nDofNod); 
        N_matrix(1, 1:2:end) = N_edge;
        N_matrix(2, 2:2:end) = N_edge;
        
        Fe = Fe + w(k) * (N_matrix' * traction) * t * J_edge;
    end
    
    dof_indices = zeros(nNodElemEdge * nDofNod, 1);
    for node_k = 1:nNodElemEdge
        node_id = edge_nodes_indices(node_k);
        dof_indices(2*node_k-1) = 2*node_id-1;
        dof_indices(2*node_k)   = 2*node_id;
    end
    
    F(dof_indices) = F(dof_indices) + Fe;
end

%% 4. BOUNDARY CONDITIONS
% MISTAKE 2 & 3 FIX: Use robust, correctly implemented boundary conditions
free = true(nDofTot, 1);
% Fix bottom edge in y-direction
bottom_edge_nodes = [1, 2, 3, 4, 5];
for node_idx = bottom_edge_nodes
    free(2 * node_idx) = false; % Correctly index the y-DOF
end
% Fix left edge in x-direction to prevent sliding and rotation
left_edge_nodes = [1, 6, 9, 14, 17];
for node_idx = left_edge_nodes
    free(2 * node_idx - 1) = false; % Correctly index the x-DOF
end

%% 5. ASSEMBLY OF GLOBAL STIFFNESS MATRIX
fprintf('2. Assembling global stiffness matrix...\n');
nEntries = nElem * (nNodElem * nDofNod)^2;
iRowSparse = zeros(nEntries, 1);
iColSparse = zeros(nEntries, 1);
valueSparseK = zeros(nEntries, 1);
iloc = 0;

for k = 1:nElem
    localNodes = mesh.nodes(mesh.elements(k, :), :);
    Kelem = getStiffnessMatrix(localNodes, 'Q8', t, E, [2 2]);

    globalDofs = zeros(1, nNodElem * nDofNod);
    for i = 1:nNodElem
        globalDofs(2*i - 1) = 2 * mesh.elements(k, i) - 1;
        globalDofs(2*i)     = 2 * mesh.elements(k, i);
    end

    for iCol = 1:(nNodElem * nDofNod)
        for iRow = 1:(nNodElem * nDofNod)
            iloc = iloc + 1;
            iRowSparse(iloc) = globalDofs(iRow);
            iColSparse(iloc) = globalDofs(iCol);
            valueSparseK(iloc) = Kelem(iRow, iCol);
        end
    end
end
KS = sparse(iRowSparse, iColSparse, valueSparseK, nDofTot, nDofTot);

%% 6. SOLVE THE SYSTEM
fprintf('3. Solving the system of equations...\n');
KS_reduced = KS(free, free);
F_reduced = F(free);
D_reduced = KS_reduced \ F_reduced;

D = zeros(nDofTot, 1);
D(free) = D_reduced;

S = stressesAtPoints(D, mesh, 'Q8', E, [2 2]);

fprintf('\n--- PATCH TEST RESULTS  ---\n');
sigma_yy_values = squeeze(S(2, :, :));

disp('Calculated sigma_yy at Gauss points for each element:');
disp(sigma_yy_values);



