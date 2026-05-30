function [p_anal, u_anal] = analyticalTerzaghi(x_query, t_query, params)
% analyticalTerzaghi  Exact series solution for 1D Terzaghi consolidation.
%   Drainage at top (x=L), impermeable at bottom (x=0).
%
%   x_query : [Nx x 1] vertical coordinates [m], x=0 at bottom
%   t_query : [Nt x 1] time values [s]
%   params  : struct with fields lambda, mu, alpha, Mbiot, kperm, muf, sigma0, L
%
%   p_anal  : [Nx x Nt]  pore pressure [Pa]
%   u_anal  : [Nx x Nt]  vertical displacement [m]  (positive = upward)
%
%   Vectorised: per time step the series is computed as a single matrix
%   product over the truncation index m (Nterms = 50).

lambda = params.lambda;
mu     = params.mu;
alpha  = params.alpha;
Mbiot  = params.Mbiot;
kperm  = params.kperm;
muf    = params.muf;
sigma0 = params.sigma0;
L      = params.L;

Moed = lambda + 2*mu;
cv   = (kperm/muf) / (1/Mbiot + alpha^2/Moed);
p0   = alpha * Mbiot * sigma0 / (Moed + alpha^2 * Mbiot);

Nterms = 50;
m      = 0:Nterms-1;
n      = 2*m + 1;                  % 1, 3, 5, ...
lam    = n * pi / (2*L);           % [1 x Nterms]
coef0  = (4/pi) * ((-1).^m ./ n);  % [1 x Nterms]  m-coefficient

x_query = x_query(:);              % [Nx x 1]
t_query = t_query(:);              % [Nt x 1]
Nx = numel(x_query);
Nt = numel(t_query);

p_anal = zeros(Nx, Nt);
u_anal = zeros(Nx, Nt);

cos_lx = cos(x_query * lam);       % [Nx x Nterms]
sin_lx = sin(x_query * lam);       % [Nx x Nterms]
inv_lam = 1 ./ lam;                % [1  x Nterms]

for it = 1:Nt
    t      = t_query(it);
    expFac = exp(-(lam.^2) * cv * t);              % [1 x Nterms]
    coeff  = coef0 .* expFac;                      % [1 x Nterms]

    p_anal(:, it) = p0 * (cos_lx * coeff(:));      % [Nx x 1]

    % u_series_m = coeff_m * sin(lam_m x) / lam_m
    u_series      = sin_lx * (coeff .* inv_lam).'; % [Nx x 1]
    u_anal(:, it) = (alpha * p0 * u_series - sigma0 * x_query) / Moed;
end
end
