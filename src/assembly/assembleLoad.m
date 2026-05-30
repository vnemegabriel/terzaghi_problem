function F = assembleLoad(mesh, eleType, sigma0, nPointsGauss)
% assembleLoad  Global mechanical load vector from uniform traction on top face.
%   Traction sigma0 [Pa] applied as compressive normal stress in -y direction
%   on all edges whose nodes all sit at y = max(y).
%
%   mesh         : struct with .nodes [nNod x 2] and .elements {nElem x 1}
%   eleType      : 'Q4' | 'Q8' | 'CST'
%   sigma0       : magnitude of applied pressure [Pa]  (positive = compression)
%   nPointsGauss : points for 1D edge integration
%   F            : [2*nNod x 1]

nNod  = size(mesh.nodes, 1);
nElem = size(mesh.elements, 1);
F     = zeros(2*nNod, 1);

ymax  = max(mesh.nodes(:,2));
tol   = 1e-10;

[w, gp] = gauss1D(nPointsGauss);

% Edge connectivity per element type (local node indices, each row is one edge)
switch eleType
    case 'Q4'
        edges = [1 2; 2 3; 3 4; 4 1];
    case 'Q8'
        edges = [1 2 5; 2 3 6; 3 4 7; 4 1 8];
    case 'CST'
        edges = [1 2; 2 3; 3 1];
    otherwise
        error('assembleLoad:unknownEleType', ...
            "Unknown eleType '%s' (expected 'Q4', 'Q8', or 'CST').", eleType)
end

% Sanity: a top face must exist in the mesh
if ~any(abs(mesh.nodes(:,2) - ymax) < tol)
    error('assembleLoad:noTopFace', ...
        'Mesh has no nodes at y = max(y) — cannot apply traction.')
end

for k = 1:nElem
    globalNodes = mesh.elements{k};
    localNodes  = mesh.nodes(globalNodes, :);

    for e = 1:size(edges, 1)
        edgeLocalIdx  = edges(e, :);
        edgeGlobalIdx = globalNodes(edgeLocalIdx);
        edgeCoords    = mesh.nodes(edgeGlobalIdx, :);

        if ~all(abs(edgeCoords(:,2) - ymax) < tol)
            continue
        end

        nNodEdge = numel(edgeLocalIdx);
        Fe = zeros(2*numel(globalNodes), 1);

        for gp_idx = 1:nPointsGauss
            zeta = gp(gp_idx);

            if nNodEdge == 2
                N_edge        = [0.5*(1-zeta), 0.5*(1+zeta)];
                dN_edge_dzeta = [-0.5, 0.5];
            else  % 3-node quadratic edge
                N_edge        = [-0.5*zeta*(1-zeta), 0.5*zeta*(1+zeta), (1-zeta^2)];
                dN_edge_dzeta = [-0.5+zeta, 0.5+zeta, -2*zeta];
            end

            dCoords_dzeta = dN_edge_dzeta * edgeCoords;
            J1D = norm(dCoords_dzeta);

            % Normal on top face: outward = +y; traction = -sigma0 * n (compression)
            t_vec = [0; -sigma0];

            for n = 1:nNodEdge
                localN  = edgeLocalIdx(n);
                dof_x   = (localN-1)*2 + 1;
                dof_y   = (localN-1)*2 + 2;
                Fe(dof_x) = Fe(dof_x) + N_edge(n) * t_vec(1) * J1D * w(gp_idx);
                Fe(dof_y) = Fe(dof_y) + N_edge(n) * t_vec(2) * J1D * w(gp_idx);
            end
        end

        % Scatter into global F
        for n = 1:numel(globalNodes)
            dof_x_g = (globalNodes(n)-1)*2 + 1;
            dof_y_g = (globalNodes(n)-1)*2 + 2;
            F(dof_x_g) = F(dof_x_g) + Fe((n-1)*2+1);
            F(dof_y_g) = F(dof_y_g) + Fe((n-1)*2+2);
        end
    end
end
end
