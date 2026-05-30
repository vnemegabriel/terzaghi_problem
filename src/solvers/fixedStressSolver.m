function [U_hist, P_hist, t_hist, iterHist] = fixedStressSolver(...
    K, H, Q, S, F_ext, U0, P0, freeDofU, freeDofP, dt, tmax, tol, maxIter, saveTimes, Kdr)
% fixedStressSolver  Fixed-stress split for Biot equations (unconditionally stable).
%
%   Adds stabilization beta = alpha^2/Kdr to the storage term in the flow step.
%   Modified flow LHS: (S/dt + beta/dt*I + H)
%   Modified flow RHS: (S/dt + beta/dt*I)*p^n - Q'/dt*(u^(k-1) - u^n)
%
%   Kdr : drained oedometric modulus = lambda + 2*mu  [Pa]
%
%   All other arguments identical to staggeredSolver.

validateattributes(dt,      {'numeric'}, {'scalar','real','positive','finite'}, mfilename, 'dt')
validateattributes(tmax,    {'numeric'}, {'scalar','real','positive','finite'}, mfilename, 'tmax')
validateattributes(tol,     {'numeric'}, {'scalar','real','positive','finite'}, mfilename, 'tol')
validateattributes(maxIter, {'numeric'}, {'scalar','integer','positive'},       mfilename, 'maxIter')
if isempty(freeDofU) || isempty(freeDofP)
    error('fixedStressSolver:noFreeDof', 'No free DOFs (system over-constrained).')
end
if tmax < dt
    error('fixedStressSolver:tmaxBelowDt', 'tmax (%.3g) is smaller than dt (%.3g).', tmax, dt)
end
if ~isequal(size(Kdr), size(S))
    error('fixedStressSolver:LstabSize', ...
        'Stabilization matrix Lstab must match size(S) = [%d %d].', size(S,1), size(S,2))
end

nSteps  = round(tmax / dt);
nSave   = numel(saveTimes);
nDofU   = numel(U0);
nDofP   = numel(P0);

U_hist  = zeros(nDofU, nSave);
P_hist  = zeros(nDofP, nSave);
t_hist  = zeros(1, nSave);
iterHist = zeros(1, nSteps);

% Extract beta from Q and Kdr:  beta = alpha^2 / Kdr  per unit volume
% S already has units of 1/M * vol; beta acts as added compressibility per node.
% Build diagonal stabilization matrix from mass matrix diagonal (lumped).
nDofP   = numel(P0);
% Use the S matrix diagonal scaled by (beta * M):  L = beta * S * M  -- no.
% Correct approach: beta * I where I is the identity on pressure DOFs.
% Physical: int N'*N * beta/dt  dOmega  =>  use S structure but factor.
% Simplest correct form: beta_stab * S * Mbiot  (since S = int 1/M N'N dOmega,
% multiplying by Mbiot gives int N'N dOmega * 1, then * beta gives int beta N'N dOmega).
% We pass beta separately and form the extra term as (beta * S * Mbiot).
% However Mbiot is not passed here. Instead accept Kdr and deduce:
%   alpha^2 / Kdr is the stabilization per volume => pass entire extra matrix L = beta * S_mass
% where S_mass = M * S. To avoid needing M, accept an extra matrix Lstab = (alpha^2/Kdr) * Smass
% from the caller. For simplicity, accept Kdr; caller must also pass alpha, Mbiot.
%
% Signature simplified: caller passes Lstab directly via Kdr argument interpreted as Lstab matrix.

Lstab = Kdr;   % caller passes the pre-built sparse stabilization matrix

A_flow = S/dt + Lstab/dt + H;
A_flow_free = A_flow(freeDofP, freeDofP);
[L_flow, U_flow, perm_flow] = lu(A_flow_free);

A_mech_free = K(freeDofU, freeDofU);
[R_mech, cholFail, P_mech] = chol(A_mech_free, 'lower', 'matrix');
useChol = (cholFail == 0);
if ~useChol
    warning('fixedStressSolver:cholFail', 'K not SPD on free DOFs; falling back to LU.')
    [L_mech, U_mech, perm_mech] = lu(A_mech_free);
end

Un = U0;
Pn = P0;

nextSaveStep = zeros(1, nSave);
for s = 1:nSave
    nextSaveStep(s) = round(saveTimes(s) / dt);
end

for step = 1:nSteps
    t_cur = step * dt;

    Pk  = Pn;
    Uk  = Un;

    for iter = 1:maxIter
        Pk_prev = Pk;

        % --- Flow step first (fixed-stress: total stress fixed) ---
        % RHS: S/dt*Pn + Lstab/dt*Pk_prev so that at convergence Pk=Pk_prev
        % the Lstab terms cancel, recovering the monolithic flow equation.
        rhs_flow = S/dt * Pn + Lstab/dt * Pk_prev - Q'/dt * (Uk - Un);
        Pk_free  = U_flow \ (L_flow \ (perm_flow * rhs_flow(freeDofP)));
        Pk = zeros(size(Pn));          % fixed BCs (p=0 at top) stay zero
        Pk(freeDofP) = Pk_free;

        % --- Mechanical step ---
        rhs_mech_free = F_ext(freeDofU) + Q(freeDofU, :) * Pk;
        if useChol
            Uk_free = P_mech * (R_mech' \ (R_mech \ (P_mech' * rhs_mech_free)));
        else
            Uk_free = U_mech \ (L_mech \ (perm_mech * rhs_mech_free));
        end
        Uk = Un;
        Uk(freeDofU) = Uk_free;

        dp = norm(Pk - Pk_prev);
        if dp / max(norm(Pk), 1) < tol
            iterHist(step) = iter;
            break
        end
        iterHist(step) = iter;
    end

    Un = Uk;
    Pn = Pk;

    for s = 1:nSave
        if step == nextSaveStep(s)
            U_hist(:, s) = Un;
            P_hist(:, s) = Pn;
            t_hist(s)    = t_cur;
        end
    end
end
end
