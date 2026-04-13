function P = eq_P_from_F(a, b, c, F)
    % Compute the 1st PK stress P from deformation gradient F and model parameters a, b, c.
    % a = dW/dI1, b = dW/dI2 at current invariants.
    % c = dW/dJ at current J.

    assert(all(size(F) == [3 3]), 'F must be 3x3');
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

    P = P_iso + P_vol; % total stress
end