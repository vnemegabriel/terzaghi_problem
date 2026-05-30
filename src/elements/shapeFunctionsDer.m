function dN = shapeFunctionsDer(locationGaussPoint, eleType)
% shapeFunctionsDer  Shape function derivatives at a single Gauss point.
%   locationGaussPoint : [1x2] [xi, eta]
%   eleType            : 'Q4' | 'Q8' | 'CST'
%   dN                 : [2 x nNodEle]  rows: d/dxi, d/deta

if numel(locationGaussPoint) < 2
    error('shapeFunctionsDer:badLocation', ...
        'locationGaussPoint must have at least 2 components [xi, eta].')
end

xi  = locationGaussPoint(1);
eta = locationGaussPoint(2);

switch eleType

    case 'Q4'
        dN1_dxi  = -0.25*(1-eta);  dN1_deta = -0.25*(1-xi);
        dN2_dxi  =  0.25*(1-eta);  dN2_deta = -0.25*(1+xi);
        dN3_dxi  =  0.25*(1+eta);  dN3_deta =  0.25*(1+xi);
        dN4_dxi  = -0.25*(1+eta);  dN4_deta =  0.25*(1-xi);
        dN = [dN1_dxi  dN2_dxi  dN3_dxi  dN4_dxi ; ...
              dN1_deta dN2_deta dN3_deta dN4_deta];

    case 'Q8'
        dN8_dxi  =  0.5*(-1)*(1-eta^2);
        dN7_dxi  =  0.5*(-2*xi)*(1+eta);
        dN6_dxi  =  0.5*(1)*(1-eta^2);
        dN5_dxi  =  0.5*(-2*xi)*(1-eta);
        dN4_dxi  =  0.25*(-1)*(1+eta) - 0.5*(dN7_dxi+dN8_dxi);
        dN3_dxi  =  0.25*(1)*(1+eta)  - 0.5*(dN6_dxi+dN7_dxi);
        dN2_dxi  =  0.25*(1)*(1-eta)  - 0.5*(dN5_dxi+dN6_dxi);
        dN1_dxi  =  0.25*(-1)*(1-eta) - 0.5*(dN5_dxi+dN8_dxi);

        dN8_deta =  0.5*(-2*eta)*(1-xi);
        dN7_deta =  0.5*(1)*(1-xi^2);
        dN6_deta =  0.5*(-2*eta)*(1+xi);
        dN5_deta =  0.5*(-1)*(1-xi^2);
        dN4_deta =  0.25*(1)*(1-xi)   - 0.5*(dN7_deta+dN8_deta);
        dN3_deta =  0.25*(1)*(1+xi)   - 0.5*(dN6_deta+dN7_deta);
        dN2_deta =  0.25*(-1)*(1+xi)  - 0.5*(dN5_deta+dN6_deta);
        dN1_deta =  0.25*(-1)*(1-xi)  - 0.5*(dN5_deta+dN8_deta);

        dN = [dN1_dxi  dN2_dxi  dN3_dxi  dN4_dxi  dN5_dxi  dN6_dxi  dN7_dxi  dN8_dxi ; ...
              dN1_deta dN2_deta dN3_deta dN4_deta dN5_deta dN6_deta dN7_deta dN8_deta];

    case 'CST'
        dN = [1  0 -1; ...
              0  1 -1];

    otherwise
        error('shapeFunctionsDer:unknownEleType', ...
            "Unknown eleType '%s' (expected 'Q4', 'Q8', or 'CST').", eleType)
end
end
