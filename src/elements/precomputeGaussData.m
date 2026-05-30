function gd = precomputeGaussData(eleType, npg)
% precomputeGaussData  Cache shape function values and derivatives at every
%   Gauss point for an element type. Reused across all elements during
%   global assembly (independent of physical coordinates).
%
%   eleType : 'Q4' | 'Q8' | 'CST'
%   npg     : number of Gauss points per direction (for Q4/Q8 only)
%
%   gd.N    : cell{nGP}    each [1 x nNodEle]
%   gd.dN   : cell{nGP}    each [2 x nNodEle]
%   gd.w    : [nGP x 1]    quadrature weights (combined wxi*weta for 2D)
%   gd.nGP  : scalar

switch eleType

    case {'Q4', 'Q8'}
        [w1, gp1] = gauss1D(npg);
        nGP = npg^2;
        gd.N  = cell(nGP, 1);
        gd.dN = cell(nGP, 1);
        gd.w  = zeros(nGP, 1);
        k = 0;
        for i = 1:npg
            for j = 1:npg
                k = k + 1;
                xi = gp1(i); eta = gp1(j);
                gd.N{k}  = shapeFunctions([xi eta], eleType);
                gd.dN{k} = shapeFunctionsDer([xi eta], eleType);
                gd.w(k)  = w1(i) * w1(j);
            end
        end
        gd.nGP = nGP;

    case 'CST'
        % Closed-form integration; single "Gauss point" at the centroid.
        gd.N   = { shapeFunctions([1/3 1/3], 'CST') };
        gd.dN  = { shapeFunctionsDer([0 0],   'CST') };
        gd.w   = 0.5;          % triangle area in reference coords
        gd.nGP = 1;

    otherwise
        error('precomputeGaussData:unknownEleType', ...
            "Unknown eleType '%s'.", eleType)
end
end
