function [cornerNodes, nodeToP, nDofP] = pressureNodeMap(mesh, eleType)
% pressureNodeMap  Corner-node list and node->pressure-DOF map for the mixed
%   u-p (Biot) formulation. Pressure is interpolated on corner nodes only
%   (Q4 basis) for both Q4 and Q8 displacement meshes, so the pressure system
%   is smaller than the displacement system on Q8 meshes. This is the single
%   source of truth used by assembly, BCs, ICs and post-processing.
%
%   eleType : 'Q4' | 'Q8' | 'CST'  (defaults to mesh.eleType if omitted)
%
%   cornerNodes : [nDofP x 1]  global indices of pressure (corner) nodes
%   nodeToP     : [nNod x 1]    map global node -> pressure DOF (0 if none)
%   nDofP       : scalar        number of pressure DOFs

if nargin < 2 || isempty(eleType)
    if isfield(mesh, 'eleType')
        eleType = mesh.eleType;
    else
        error('pressureNodeMap:missingEleType', ...
            'eleType not supplied and not present in mesh struct.')
    end
end

% Corner nodes are the leading local nodes of each element's connectivity:
%   Q4/Q8 -> first 4 (the geometric corners), CST -> all 3.
switch eleType
    case {'Q4', 'Q8'}
        nLocCorner = 4;
    case 'CST'
        nLocCorner = 3;
    otherwise
        error('pressureNodeMap:unknownEleType', ...
            "Unknown eleType '%s'.", eleType)
end

nNod     = size(mesh.nodes, 1);
isCorner = false(nNod, 1);
for k = 1:numel(mesh.elements)
    en = mesh.elements{k};
    isCorner(en(1:nLocCorner)) = true;
end

cornerNodes        = find(isCorner);
nDofP              = numel(cornerNodes);
nodeToP            = zeros(nNod, 1);
nodeToP(cornerNodes) = 1:nDofP;
end
