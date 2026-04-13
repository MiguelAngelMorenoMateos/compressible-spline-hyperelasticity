function rvec = loss_function(x, model, prot, data)
    % Residual function for lsqlin:
    %   minimizes ||rvec||^2

    % Count total residuals across all modes
    Ntot = 0;
    for i = 1:numel(data.modes)
        mode = data.modes{i};
        Ntot = Ntot + numel(data.(mode).P11);
    end

    rvec = zeros(Ntot, 1);
    rvecP22 = [];

    % return nan if any parameter is nan or inf
    if any(isnan(x)) || any(isinf(x))
        rvec(:) = nan;
        return;
    end

    model = model.unpack(x, model); % update splines

    row = 0;
    for i = 1:numel(data.modes) % UT, UC, SS
        mode = data.modes{i};
        D = data.(mode);
        p = prot.(mode);

        Nm = numel(D.P11);
        Pmax = max(abs(D.P11));
        scale1 = sqrt(D.w/Nm/Pmax);          % For experiment-model fitting in the lsq problems.
        scale2 = 100*sqrt(D.w/Nm/Pmax);      % For P22=0 BC in the lsq problems.

        for k = 1:Nm % Loop over data points in this mode (all strain-stress pairs of experimental data)
            row = row + 1;
            lam = D.lambda(k); 
            F   = p.Fbar(lam); 
            [a,b,c] = model.eval_ab(model, F);
            P_pred = eq_P_from_F(a,b,c,F);

            if strcmpi(mode,'SS') % For Shear test (SS): residual is r = P12 - data
                rvec(row) =  scale1 * (P_pred(1,2) - D.P11(k)); % Difference between analytical and experimental stress component. (Note, D.P11 is the name of the variable in the data struct, but it actually contains the P12 data for the shear test. This is a bit confusing but we can live with it.)
            else % For Uniaxial test (UT & UC): residual is r = P11 - data
                rvec(row) =  scale1 * (P_pred(1,1) - D.P11(k)); % Difference between analytical and experimental stress component.
                rvecP22(end+1) = scale2 * P_pred(2,2); % residual to enforce zero lateral stress for uniaxial tests.
            end
        end
    end
    rvec = [rvec; rvecP22']; % Append P22 residuals for uniaxial tests to the end of the residual vector.
end