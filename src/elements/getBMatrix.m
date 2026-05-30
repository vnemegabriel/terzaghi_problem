function [B, detJ] = getBMatrix(locationGaussPoint, nodesElem, eleType)
% getBMatrix  Strain-displacement matrix B and Jacobian determinant.
%   locationGaussPoint : [1x2] [xi, eta]
%   nodesElem          : [nNodEle x 2] global coordinates
%   eleType            : 'Q4' | 'Q8' | 'CST'
%   B                  : [3 x 2*nNodEle]
%   detJ               : scalar

nNodEle = size(nodesElem, 1);
nDofEle = 2 * nNodEle;

dN_iso = shapeFunctionsDer(locationGaussPoint, eleType);
Jacobian = dN_iso * nodesElem;
detJ = det(Jacobian);
dNxy = Jacobian \ dN_iso;   % [2 x nNodEle]  d/dx, d/dy

B = zeros(3, nDofEle);
B(1, 1:2:end) = dNxy(1,:);
B(2, 2:2:end) = dNxy(2,:);
B(3, 1:2:end) = dNxy(2,:);
B(3, 2:2:end) = dNxy(1,:);
end
