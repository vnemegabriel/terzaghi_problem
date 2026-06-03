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
%   Mixed u-p (Biot) formulation: displacement and pressure are NOT
%   interpolated with the same shape functions. The strain-displacement
%   matrix B and the geometry Jacobian use the full displacement basis
%   (gd.N / gd.dN), while the pressure field uses the corner-node basis
%   (gd.Np / gd.dNp) — Q4 for both Q4 and Q8 elements (Taylor-Hood). This
%   satisfies the inf-sup (LBB) stability condition for the coupled system.
%
%   nP = gd.nP = number of pressure (corner) nodes per element
%               (4 for Q4/Q8, 3 for CST).
%
%   He : [nP x nP]          permeability matrix
%   Qe : [2*nNodEle x nP]   coupling matrix
%   Se : [nP x nP]          storage matrix

if nargin < 8 || isempty(gd)
    gd = precomputeGaussData(eleType, nPointsGauss);
end

nDofNod = 2;
nNodEle = size(nodesElem, 1);
nDofEle = nDofNod * nNodEle;
nP      = gd.nP;                   % pressure (corner) nodes

m_vec = [1; 1; 0];                 % Voigt volumetric vector
kf    = kperm / muf;

He = zeros(nP);
Qe = zeros(nDofEle, nP);
Se = zeros(nP);

switch eleType

    case {'Q4', 'Q8'}
        for k = 1:gd.nGP
            dN_iso  = gd.dN{k};    % displacement / geometry basis
            Np      = gd.Np{k};    % pressure shape functions (corner)  [1 x nP]
            dNp_iso = gd.dNp{k};   % pressure shape func derivatives     [2 x nP]

            % Geometry mapped by the displacement (full) basis
            Jacobian = dN_iso * nodesElem;
            assert(det(Jacobian) > 0)

            dNxy  = Jacobian \ dN_iso;     % displacement gradients (for B)
            dNpxy = Jacobian \ dNp_iso;    % pressure gradients

            B = zeros(3, nDofEle);
            B(1, 1:2:end) = dNxy(1,:);
            B(2, 2:2:end) = dNxy(2,:);
            B(3, 1:2:end) = dNxy(2,:);
            B(3, 2:2:end) = dNxy(1,:);

            fac = gd.w(k) * det(Jacobian);
            He = He + kf  * (dNpxy' * dNpxy)    * fac;
            Qe = Qe + alpha * B' * m_vec * Np    * fac;
            Se = Se + (1/Mbiot) * (Np' * Np)     * fac;
        end

    case 'CST'
        dN_iso  = gd.dN{1};
        Np      = gd.Np{1};        % corner basis == displacement basis (linear)
        dNp_iso = gd.dNp{1};

        Jacobian = dN_iso * nodesElem;
        detJ = det(Jacobian);
        Area = detJ / 2;

        dNxy  = Jacobian \ dN_iso;
        dNpxy = Jacobian \ dNp_iso;

        B = zeros(3, nDofEle);
        B(1, 1:2:end) = dNxy(1,:);
        B(2, 2:2:end) = dNxy(2,:);
        B(3, 1:2:end) = dNxy(2,:);
        B(3, 2:2:end) = dNxy(1,:);

        He = kf  * (dNpxy' * dNpxy)    * Area;
        Qe = alpha * B' * m_vec * Np    * Area;
        Se = (1/Mbiot) * (Np' * Np)     * Area;

    otherwise
        error('assemblePoroelastic:unknownEleType', ...
            "Unknown eleType '%s'.", eleType)
end
end
