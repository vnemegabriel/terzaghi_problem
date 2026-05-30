function [w, gp] = gauss1D(n)
% gauss1D  1-D Gauss-Legendre quadrature weights and points on [-1, 1].
%   n   : number of points (1, 2, or 3)
%   w   : column vector of weights
%   gp  : column vector of integration points

switch n
    case 1
        w  = 2;
        gp = 0;
    case 2
        w  = [1; 1];
        a  = sqrt(3)/3;
        gp = [-a; a];
    case 3
        w  = [5; 8; 5] / 9;
        a  = sqrt(3/5);
        gp = [-a; 0; a];
    otherwise
        error('gauss1D:invalidN', 'Only n = 1, 2, or 3 supported.')
end
end
