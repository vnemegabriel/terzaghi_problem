function [He, Qe, Se] = assemblePoroelastic(nodesElem, eleType, kperm, muf, alpha, Mbiot, nPointsGauss, gd)
% assemblePoroelastic  Element poroelastic matrices in one Gauss pass.
%   nodesElem    : [nNodEle x 2] global coordinates
%   eleType      : 'Q4' | 'Q8' | 'CST'
%   kperm        : intrinsic permeability [m^2]
%   muf          : fluid viscosity [Pa*s]
%   alpha        : Biot coefficient [-]
%   Mbiot        : Biot modulus [Pa]
%   nPointsGauss : scalar
%   gd           : (optional) precomputed shape data from precomputeGaussData
%
%   He : [nNodEle x nNodEle]    permeability matrix
%   Qe : [2*nNodEle x nNodEle]  coupling matrix
%   Se : [nNodEle x nNodEle]    storage matrix

if nargin < 8 || isempty(gd)
    gd = precomputeGaussData(eleType, nPointsGauss);
end

nDofNod = 2;
nNodEle = size(nodesElem, 1);
nDofEle = nDofNod * nNodEle;

m_vec = [1; 1; 0];                 % Voigt volumetric vector
kf    = kperm / muf;

He = zeros(nNodEle);
Qe = zeros(nDofEle, nNodEle);
Se = zeros(nNodEle);

switch eleType

    case {'Q4', 'Q8'}
        for k = 1:gd.nGP
            N_vec  = gd.N{k};
            dN_iso = gd.dN{k};

            Jacobian = dN_iso * nodesElem;
            assert(det(Jacobian) > 0)

            dNxy = Jacobian \ dN_iso;

            B = zeros(3, nDofEle);
            B(1, 1:2:end) = dNxy(1,:);
            B(2, 2:2:end) = dNxy(2,:);
            B(3, 1:2:end) = dNxy(2,:);
            B(3, 2:2:end) = dNxy(1,:);

            fac = gd.w(k) * det(Jacobian);
            He = He + kf  * (dNxy' * dNxy)         * fac;
            Qe = Qe + alpha * B' * m_vec * N_vec    * fac;
            Se = Se + (1/Mbiot) * (N_vec' * N_vec)  * fac;
        end

    case 'CST'
        N_vec  = gd.N{1};
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

        He = kf  * (dNxy' * dNxy)        * Area;
        Qe = alpha * B' * m_vec * N_vec   * Area;
        Se = (1/Mbiot) * (N_vec' * N_vec) * Area;

    otherwise
        error('assemblePoroelastic:unknownEleType', ...
            "Unknown eleType '%s'.", eleType)
end
end
