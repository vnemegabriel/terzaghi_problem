function [K, H, Q, S] = buildGlobalMatrices(mesh, eleType, params, npg)
% buildGlobalMatrices  Assemble global K, H, Q, S using triplet sparse format.

    if nargin < 2 || isempty(eleType)
        if isfield(mesh, 'eleType')
            eleType = mesh.eleType;
        else
            error('buildGlobalMatrices:missingEleType', ...
                'eleType not supplied and not present in mesh struct.')
        end
    end

    nNod  = size(mesh.nodes, 1);
    nElem = size(mesh.elements, 1);

    % Plane-strain constitutive matrix (default mode)
    C = getConstitutiveMatrix(params.lambda, params.mu, 'plane_strain');

    % Mixed u-p: displacement on all nodes, pressure on corner nodes only.
    [~, nodeToP, nDofP] = pressureNodeMap(mesh, eleType);
    nDofU = 2*nNod;

    % Estimate element size from first element (assumes uniform mesh)
    nNodEleMax = max(cellfun(@numel, mesh.elements));
    nDofEleMax = 2 * nNodEleMax;

    % Pre-allocate triplet storage
    nK = nElem * nDofEleMax^2;
    nH = nElem * nNodEleMax^2;
    nQ = nElem * nDofEleMax * nNodEleMax;
    nS = nH;

    iK = zeros(nK,1);  jK = zeros(nK,1);  vK = zeros(nK,1);
    iH = zeros(nH,1);  jH = zeros(nH,1);  vH = zeros(nH,1);
    iQ = zeros(nQ,1);  jQ = zeros(nQ,1);  vQ = zeros(nQ,1);
    iS = zeros(nS,1);  jS = zeros(nS,1);  vS = zeros(nS,1);

    pK = 0; pH = 0; pQ = 0; pS = 0;

    % Pre-compute shape data once for the chosen eleType / npg
    gd = precomputeGaussData(eleType, npg);

    for k = 1:nElem
        globalNodes = mesh.elements{k};
        localNodes  = mesh.nodes(globalNodes, :);
        nNodEle     = numel(globalNodes);

        Ke = assembleMechanical(localNodes, eleType, C, npg, gd);
        [He, Qe, Se] = assemblePoroelastic(localNodes, eleType, ...
            params.kperm, params.muf, params.alpha, params.Mbiot, npg, gd);

        % Mechanical DOF map (all element nodes)
        dofU = zeros(1, 2*nNodEle);
        dofU(1:2:end) = 2*globalNodes - 1;
        dofU(2:2:end) = 2*globalNodes;

        % Pressure DOF map (corner nodes only -> pressure DOF space).
        % He/Se are [nP x nP], Qe is [2*nNodEle x nP] with nP corner nodes.
        nP   = size(He, 1);
        dofP = nodeToP(globalNodes(1:nP));
        dofP = dofP(:)';

        % Scatter Ke
        [II, JJ] = ndgrid(dofU, dofU);
        m = numel(II);
        iK(pK+1:pK+m) = II(:);  jK(pK+1:pK+m) = JJ(:);  vK(pK+1:pK+m) = Ke(:);
        pK = pK + m;

        % Scatter He
        [II, JJ] = ndgrid(dofP, dofP);
        m = numel(II);
        iH(pH+1:pH+m) = II(:);  jH(pH+1:pH+m) = JJ(:);  vH(pH+1:pH+m) = He(:);
        iS(pS+1:pS+m) = II(:);  jS(pS+1:pS+m) = JJ(:);  vS(pS+1:pS+m) = Se(:);
        pH = pH + m;  pS = pS + m;

        % Scatter Qe (rectangular)
        [II, JJ] = ndgrid(dofU, dofP);
        m = numel(II);
        iQ(pQ+1:pQ+m) = II(:);  jQ(pQ+1:pQ+m) = JJ(:);  vQ(pQ+1:pQ+m) = Qe(:);
        pQ = pQ + m;
    end

    % Trim trailing zeros (in case some elements were smaller than max)
    iK = iK(1:pK); jK = jK(1:pK); vK = vK(1:pK);
    iH = iH(1:pH); jH = jH(1:pH); vH = vH(1:pH);
    iQ = iQ(1:pQ); jQ = jQ(1:pQ); vQ = vQ(1:pQ);
    iS = iS(1:pS); jS = jS(1:pS); vS = vS(1:pS);

    K = sparse(iK, jK, vK, nDofU, nDofU);
    H = sparse(iH, jH, vH, nDofP, nDofP);
    Q = sparse(iQ, jQ, vQ, nDofU, nDofP);
    S = sparse(iS, jS, vS, nDofP, nDofP);
end
