function gd = precomputeGaussData(eleType, npg)
% precomputeGaussData  Cache shape function values and derivatives at every
%   Gauss point for an element type. Reused across all elements during
%   global assembly (independent of physical coordinates).
%
%   eleType : 'Q4' | 'Q8' | 'CST'
%   npg     : number of Gauss points per direction (for Q4/Q8 only)
%
%   Two bases are cached at the same Gauss points:
%     * displacement / geometry basis (full eleType) — used for K, B, and the
%       isoparametric geometry map (Jacobian);
%     * pressure basis (corner nodes only) — used for the poroelastic H, Q, S.
%   These are kept SEPARATE on purpose: the mixed u-p (Biot) formulation must
%   not interpolate displacement and pressure with the same shape functions
%   (Taylor-Hood / inf-sup requirement). For Q8 the pressure basis is the
%   corner (Q4) basis; for Q4 it coincides with the displacement basis.
%
%   gd.N    : cell{nGP}    each [1 x nNodEle]  displacement shape functions
%   gd.dN   : cell{nGP}    each [2 x nNodEle]  displacement shape func derivs
%   gd.Np   : cell{nGP}    each [1 x nP]       pressure (corner) shape functions
%   gd.dNp  : cell{nGP}    each [2 x nP]       pressure (corner) shape func derivs
%   gd.w    : [nGP x 1]    quadrature weights (combined wxi*weta for 2D)
%   gd.nGP  : scalar       number of Gauss points
%   gd.nP   : scalar       number of pressure (corner) nodes per element

switch eleType

    case {'Q4', 'Q8'}
        % Pressure is interpolated on the corner nodes only (bilinear Q4
        % basis) for BOTH Q4 and Q8 displacement elements.
        pEleType = 'Q4';

        [w1, gp1] = gauss1D(npg);
        nGP = npg^2;
        gd.N   = cell(nGP, 1);
        gd.dN  = cell(nGP, 1);
        gd.Np  = cell(nGP, 1);
        gd.dNp = cell(nGP, 1);
        gd.w   = zeros(nGP, 1);
        k = 0;
        for i = 1:npg
            for j = 1:npg
                k = k + 1;
                xi = gp1(i); eta = gp1(j);
                gd.N{k}   = shapeFunctions([xi eta], eleType);
                gd.dN{k}  = shapeFunctionsDer([xi eta], eleType);
                gd.Np{k}  = shapeFunctions([xi eta], pEleType);
                gd.dNp{k} = shapeFunctionsDer([xi eta], pEleType);
                gd.w(k)   = w1(i) * w1(j);
            end
        end
        gd.nGP = nGP;
        gd.nP  = 4;            % corner nodes (Q4 pressure basis)

    case 'CST'
        % Closed-form integration; single "Gauss point" at the centroid.
        % The linear CST basis is already corner-only, so pressure and
        % displacement bases coincide.
        gd.N   = { shapeFunctions([1/3 1/3], 'CST') };
        gd.dN  = { shapeFunctionsDer([0 0],   'CST') };
        gd.Np  = gd.N;
        gd.dNp = gd.dN;
        gd.w   = 0.5;          % triangle area in reference coords
        gd.nGP = 1;
        gd.nP  = 3;

    otherwise
        error('precomputeGaussData:unknownEleType', ...
            "Unknown eleType '%s'.", eleType)
end
end
