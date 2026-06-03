function [U0, P0] = undrainedIC(K, Q, F_ext, freeDofU, freeDofP, p0)
    % Undrained: p = p0 at interior (free) DOFs, p=0 at drained boundary.
    % Mechanical equilibrium solved with this pressure field.
    % Sizes are inferred from the coupling matrix Q [nDofU x nDofP], so this
    % works for both equal-order (Q4) and corner-pressure (Q8) systems.
    nDofU = size(Q, 1);
    nDofP = size(Q, 2);
    P0 = zeros(nDofP, 1);
    P0(freeDofP) = p0;          % p=0 enforced at fixed (top) pressure DOFs
    rhs = F_ext + Q * P0;
    U0  = zeros(nDofU, 1);
    U0(freeDofU) = K(freeDofU, freeDofU) \ rhs(freeDofU);
end

