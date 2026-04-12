function dP_dx = eq_dPdx_from_F_with_IFT(model, da_dx, db_dx, dc_dx, F, lam2_free)
% Total derivative d vec(P) / d x when lateral stretch lam is defined implicitly by P22=0.
% Assumes F(2,2)=F(3,3)=lam at the converged state.

    % --- 1) explicit partial at fixed F ---
    dPdx_fixed = eq_dPdx_from_F(da_dx, db_dx, dc_dx, F);   % 9 x nParams

    % Indices in MATLAB column-major P(:):
    idx11 = 1;
    idx22 = 5;

    % --- 2) compute dP/dlam (finite diff in lam only) ---
    lam = F(2,2);
    h = 1e-6 * max(1, abs(lam));

    Fp = F; Fp(2,2) = lam + h; Fp(3,3) = lam + h;
    Fm = F; Fm(2,2) = lam - h; Fm(3,3) = lam - h;

    [ap,bp,cp] = model.eval_ab(model, Fp);
    Pp = eq_P_from_F(ap,bp,cp,Fp);

    [am,bm,cm] = model.eval_ab(model, Fm);
    Pm = eq_P_from_F(am,bm,cm,Fm);
    dP_dlam = (Pp(:) - Pm(:)) / (2*h);   % 9 x 1

    % --- 3) implicit function theorem for dlam/dx ---
    % g(lam,x) = P22(lam,x) = 0
    % dg/dx = dP22/dx at fixed F  (but with lam held -> exactly idx22 row)
    dg_dx = dPdx_fixed(idx22, :);        % 1 x nParams

    dg_dlam = dP_dlam(idx22);            % scalar

    if abs(dg_dlam) < 1e-14
        warning('dg_dlam ~ 0; implicit correction ill-conditioned. Returning fixed-F derivative.');
        dP_dx = dPdx_fixed;
        return;
    end

    dlam_dx = - dg_dx / dg_dlam;         % 1 x nParams

    % --- 4) chain rule: dP/dx = dP/dx|F + (dP/dlam) * (dlam/dx) ---
    if lam2_free
        dP_dx = dPdx_fixed + dP_dlam * dlam_dx;   % (9x1)*(1xn) => 9 x nParams  % NOTE: second term only necessary if lam2 is free variable!
    else
        dP_dx = dPdx_fixed; % If lam2 is fixed, then no implicit dependence of lam on parameters, so just return fixed-F derivative.
    end
    
end
