function Ni = shapeFunctions(locationGaussPoint, eleType)
% shapeFunctions  Shape function values at a single Gauss point.
%   locationGaussPoint : [1x2] [xi, eta]
%   eleType            : 'Q4' | 'Q8' | 'CST'
%   Ni                 : [1 x nNodEle]

if numel(locationGaussPoint) < 2
    error('shapeFunctions:badLocation', ...
        'locationGaussPoint must have at least 2 components [xi, eta].')
end

xi  = locationGaussPoint(1);
eta = locationGaussPoint(2);

switch eleType

    case 'Q4'
        N1 = 0.25*(1-xi)*(1-eta);
        N2 = 0.25*(1+xi)*(1-eta);
        N3 = 0.25*(1+xi)*(1+eta);
        N4 = 0.25*(1-xi)*(1+eta);
        Ni = [N1 N2 N3 N4];

    case 'Q8'
        N8 = 0.5*(1-xi)*(1-eta^2);
        N7 = 0.5*(1-xi^2)*(1+eta);
        N6 = 0.5*(1+xi)*(1-eta^2);
        N5 = 0.5*(1-xi^2)*(1-eta);
        N4 = 0.25*(1-xi)*(1+eta) - 0.5*(N7+N8);
        N3 = 0.25*(1+xi)*(1+eta) - 0.5*(N6+N7);
        N2 = 0.25*(1+xi)*(1-eta) - 0.5*(N5+N6);
        N1 = 0.25*(1-xi)*(1-eta) - 0.5*(N8+N5);
        Ni = [N1 N2 N3 N4 N5 N6 N7 N8];

    case 'CST'
        N1 = xi;
        N2 = eta;
        N3 = 1 - xi - eta;
        Ni = [N1 N2 N3];

    otherwise
        error('shapeFunctions:unknownEleType', ...
            "Unknown eleType '%s' (expected 'Q4', 'Q8', or 'CST').", eleType)
end
end
