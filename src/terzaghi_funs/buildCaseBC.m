function [freeDofU, freeDofP] = buildCaseBC(mesh, W, L)
    nNod = size(mesh.nodes, 1);
    tol  = 1e-10;

    fixedU = false(nNod, 2);   % [ux, uy]
    for n = 1:nNod
        x = mesh.nodes(n,1);
        y = mesh.nodes(n,2);
        if abs(x) < tol || abs(x - W) < tol    % left / right: fix ux
            fixedU(n, 1) = true;
        end
        if abs(y) < tol                          % bottom: fix uy
            fixedU(n, 2) = true;
        end
    end

    allDofU  = 1:2*nNod;
    fixedVec = reshape(fixedU', [], 1);
    freeDofU = allDofU(~fixedVec);

    % Pressure lives on corner nodes only (mixed u-p formulation), so the
    % pressure BC is expressed in pressure-DOF space, not node space.
    [cornerNodes, ~, nDofP] = pressureNodeMap(mesh);
    fixedP = abs(mesh.nodes(cornerNodes,2) - L) < tol;   % top: p=0
    allDofP  = 1:nDofP;
    freeDofP = allDofP(~fixedP);
end