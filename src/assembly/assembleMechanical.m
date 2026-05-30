function Ke = assembleMechanical(nodesElem, eleType, C, nPointsGauss, gd)
% assembleMechanical  Element mechanical stiffness matrix Ke = int B'*C*B dOmega.
%   nodesElem    : [nNodEle x 2] global coordinates
%   eleType      : 'Q4' | 'Q8' | 'CST'
%   C            : [3x3] constitutive matrix (plane strain)
%   nPointsGauss : scalar, number of Gauss points per direction
%   gd           : (optional) precomputed shape data from precomputeGaussData.
%                  If omitted, computed in-place.
%   Ke           : [2*nNodEle x 2*nNodEle]

if nargin < 5 || isempty(gd)
    gd = precomputeGaussData(eleType, nPointsGauss);
end

nDofNod = 2;
nNodEle = size(nodesElem, 1);
nDofEle = nDofNod * nNodEle;

switch eleType

    case {'Q4', 'Q8'}
        Ke = zeros(nDofEle);

        for k = 1:gd.nGP
            dN_iso = gd.dN{k};

            Jacobian = dN_iso * nodesElem;
            assert(det(Jacobian) > 0)

            dNxy = Jacobian \ dN_iso;

            % Strain-displacement B matrix
            B = zeros(3, nDofEle);
            B(1, 1:2:end) = dNxy(1,:);
            B(2, 2:2:end) = dNxy(2,:);
            B(3, 1:2:end) = dNxy(2,:);
            B(3, 2:2:end) = dNxy(1,:);

            Ke = Ke + gd.w(k) * B'*C*B * det(Jacobian);
        end

    case 'CST'
        dN_iso = gd.dN{1};

        Jacobian = dN_iso * nodesElem;
        detJ = det(Jacobian);
        Area = detJ / 2;

        dNxy = Jacobian \ dN_iso;

        B = zeros(3, nDofEle);
        B(1, 1:2:end) = dNxy(1,:);
        B(2, 2:2:end) = dNxy(2,:);
        B(3, 1:2:end) = dNxy(2,:);
        B(3, 2:2:end) = dNxy(1,:);

        Ke = B' * C * B * Area;

    otherwise
        error('assembleMechanical:unknownEleType', ...
            "Unknown eleType '%s'.", eleType)
end
end
