function [U_hist, P_hist, t_hist, iterHist, pNormHist] = staggeredSolver(...
    K, H, Q, S, F_ext, U0, P0, freeDofU, freeDofP, dt, tmax, tol, maxIter, saveTimes)
% staggeredSolver  Undrained-split staggered time integration for Biot equations.
%
%   System at each step:
%     Mechanical : K * u^(k)  = F_ext + Q * p^(k-1)
%     Flow       : (S/dt + H) * p^(k) = S/dt * p^n - Q'/dt * (u^(k) - u^n)
%
%   Convergence : ||p^(k) - p^(k-1)||_2 / max(||p^(k)||_2, 1) < tol
%
%   K, H, Q, S   : global matrices (sparse)
%   F_ext        : [nDofU x 1]  constant external force (traction)
%   U0, P0       : initial conditions
%   freeDofU/P   : logical or index vectors for free DOFs
%   dt           : time step [s]
%   tmax         : end time [s]
%   tol          : convergence tolerance
%   maxIter      : max inner iterations per step
%   saveTimes    : [Ns x 1] times at which to store results (snapped to nearest step)
%
%   U_hist       : [nDofU x Ns]
%   P_hist       : [nDofP x Ns]
%   t_hist       : [1 x Ns]  actual saved times
%   iterHist     : [1 x nSteps]  iterations per step
%   pNormHist    : {nSteps x 1}  norm of p per inner iteration (last step of divergence study)

validateattributes(dt,      {'numeric'}, {'scalar','real','positive','finite'}, mfilename, 'dt')
validateattributes(tmax,    {'numeric'}, {'scalar','real','positive','finite'}, mfilename, 'tmax')
validateattributes(tol,     {'numeric'}, {'scalar','real','positive','finite'}, mfilename, 'tol')
validateattributes(maxIter, {'numeric'}, {'scalar','integer','positive'},       mfilename, 'maxIter')
if isempty(freeDofU)
    error('staggeredSolver:noFreeMechDof', 'No free mechanical DOFs (system over-constrained).')
end
if isempty(freeDofP)
    error('staggeredSolver:noFreePressureDof', 'No free pressure DOFs (system over-constrained).')
end
if tmax < dt
    error('staggeredSolver:tmaxBelowDt', 'tmax (%.3g) is smaller than dt (%.3g).', tmax, dt)
end

nSteps  = round(tmax / dt);
nSave   = numel(saveTimes);
nDofU   = numel(U0);
nDofP   = numel(P0);

U_hist  = zeros(nDofU, nSave);
P_hist  = zeros(nDofP, nSave);
t_hist  = zeros(1, nSave);
iterHist   = zeros(1, nSteps);
pNormHist  = cell(nSteps, 1);

Un = U0;
Pn = P0;

% Precompute LHS for flow step (constant over time)
A_flow = S/dt + H;
A_flow_free = A_flow(freeDofP, freeDofP);
[L_flow, U_flow, perm_flow] = lu(A_flow_free);

% Precompute LHS for mechanical step.
% K is SPD on free DOFs -> Cholesky is 2x faster than LU and half the memory.
A_mech_free = K(freeDofU, freeDofU);
[R_mech, cholFail, P_mech] = chol(A_mech_free, 'lower', 'matrix');
useChol = (cholFail == 0);
if ~useChol
    % Fallback to LU if Cholesky declines (e.g. non-SPD due to BC pattern)
    warning('staggeredSolver:cholFail', 'K not SPD on free DOFs; falling back to LU.')
    [L_mech, U_mech, perm_mech] = lu(A_mech_free);
end

saveIdx = 0;
nextSaveStep = zeros(1, nSave);
for s = 1:nSave
    nextSaveStep(s) = round(saveTimes(s) / dt);
end

for step = 1:nSteps
    t_cur = step * dt;

    Pk  = Pn;  % initial guess: pressure from previous step
    Uk  = Un;
    pNorms_step = zeros(maxIter, 1);
    converged   = false;

    for iter = 1:maxIter
        Pk_prev = Pk;

        % --- Mechanical step: fix Pk ---
        rhs_mech = F_ext + Q * Pk;
        b        = rhs_mech(freeDofU);
        if useChol
            % Solve P*K*P' * (P*x) = P*b   with R_mech = chol(P*K*P','lower')
            Uk_free = P_mech * (R_mech' \ (R_mech \ (P_mech' * b)));
        else
            Uk_free = U_mech \ (L_mech \ (perm_mech * b));
        end
        Uk = Un;
        Uk(freeDofU) = Uk_free;

        % --- Flow step: fix Uk ---
        rhs_flow = S/dt * Pn - Q'/dt * (Uk - Un);
        Pk_free  = U_flow \ (L_flow \ (perm_flow * rhs_flow(freeDofP)));
        Pk = zeros(size(Pn));          % fixed BCs (p=0 at top) stay zero
        Pk(freeDofP) = Pk_free;

        dp = norm(Pk - Pk_prev);
        pNorms_step(iter) = norm(Pk);
        if dp / max(norm(Pk), 1) < tol
            converged = true;
            iterHist(step) = iter;
            break
        end

        if ~isfinite(dp)
            iterHist(step) = iter;
            pNormHist{step} = pNorms_step(1:iter);
            fprintf('  Step %d: DIVERGED at iter %d (t=%.4g s)\n', step, iter, t_cur);
            break
        end
    end

    if ~converged && isfinite(dp)
        iterHist(step) = maxIter;
        fprintf('  Step %d: did not converge in %d iter (t=%.4g s)\n', step, maxIter, t_cur);
    end
    pNormHist{step} = pNorms_step(1:iterHist(step));

    Un = Uk;
    Pn = Pk;

    % Save if requested
    for s = 1:nSave
        if step == nextSaveStep(s)
            saveIdx = saveIdx + 1;
            U_hist(:, s) = Un;
            P_hist(:, s) = Pn;
            t_hist(s)    = t_cur;
        end
    end

    if ~isfinite(norm(Pn))
        fprintf('Stopping: solution blew up at step %d (t=%.4g s).\n', step, t_cur);
        U_hist = U_hist(:, 1:saveIdx);
        P_hist = P_hist(:, 1:saveIdx);
        t_hist = t_hist(1:saveIdx);
        return
    end
end
end
