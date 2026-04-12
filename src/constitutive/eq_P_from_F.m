function P = eq_P_from_F(a, b, c, F)
% Compute P for incompressible isotropic hyperelasticity given
% a = dW/dI1, b = dW/dI2 at current invariants.
%
% This is formulation-agnostic: caller provides a,b.

    assert(all(size(F) == [3 3]), 'F must be 3x3');
    %assert(abs(det(F) - 1) < 1e-6, 'F must be volume-preserving (det(F)=1)');

    % fzero local solver for lam2 or lam3
    %fun = @(lam) eq_P33_from_lam(lam, a, b, F);
    %lam3 = fzero(fun, 1.0);

    % update F with solved lam3
    %F(3,3) = lam3; F(2,2) = 1/lam3; % enforce incompressibility

    C = F.' * F;
    Cbar = C * det(F)^(-2/3); % isochoric part of C
    I1 = trace(C);
    I1bar = I1 * det(F)^(-2/3);
    I  = eye(3);

    % K = (I1*I - C)
    K = I1bar * I - Cbar;

    % 2nd PK (isochoric part)
    S_bar = 2 * (a * I + b * K);

    Cinv = inv(C);
    trace_term = sum(sum(C .* S_bar));        % C : S_bar
    S_iso = det(F)^(-2/3) * (S_bar - (1/3) * trace_term * Cinv);  % Eq. 5, Steinmann et al. 2012?

    % S_vol = c * dJdC
    S_vol = c * det(F) * Cinv; 

    % 1st PK
    P_iso = F * S_iso;
    P_vol = F * S_vol;

    % pressure elimination with P33=0 => p = F33*Piso33
    % p = F(3,3) * P_iso(3,3);
    % FinvT = inv(F).';

    P = P_iso + P_vol; % total stress
end