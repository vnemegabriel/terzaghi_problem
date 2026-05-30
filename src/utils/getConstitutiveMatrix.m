function C = getConstitutiveMatrix(lambda, mu, mode)
% getConstitutiveMatrix  Linear-elastic 2D constitutive matrix C [3x3].
%
%   lambda, mu : Lamé parameters [Pa]
%   mode       : 'plane_strain' (default) | 'plane_stress'
%   C          : Voigt stiffness s.t.  sigma_v = C * eps_v
%                with eps_v = [eps_xx; eps_yy; 2 eps_xy]

if nargin < 3 || isempty(mode)
    mode = 'plane_strain';
end

switch lower(mode)
    case 'plane_strain'
        C = [lambda + 2*mu,  lambda,         0 ; ...
             lambda,         lambda + 2*mu,  0 ; ...
             0,              0,              mu];

    case 'plane_stress'
        % Convert (lambda, mu) -> (E, nu) then build plane-stress C
        E  = mu * (3*lambda + 2*mu) / (lambda + mu);
        nu = lambda / (2 * (lambda + mu));
        f  = E / (1 - nu^2);
        C  = f * [1   nu  0          ; ...
                  nu  1   0          ; ...
                  0   0   (1 - nu)/2];

    otherwise
        error('getConstitutiveMatrix:unknownMode', ...
            "Unknown mode '%s' (expected 'plane_strain' or 'plane_stress').", mode)
end
end
