function [U0, P0] = undrainedIC(K, Q, F_ext, freeDofU, freeDofP, p0, nNod)
    % Undrained: p = p0 at interior (free) DOFs, p=0 at drained boundary.
    % Mechanical equilibrium solved with this pressure field.
    P0 = zeros(nNod, 1);
    P0(freeDofP) = p0;          % p=0 enforced at fixed (top) nodes
    rhs = F_ext + Q * P0;
    U0  = zeros(2*nNod, 1);
    U0(freeDofU) = K(freeDofU, freeDofU) \ rhs(freeDofU);
end

