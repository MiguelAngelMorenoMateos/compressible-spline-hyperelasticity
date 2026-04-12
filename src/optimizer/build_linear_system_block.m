function [A, y, blockIdx] = build_linear_system_block(model, prot, data, blockNames)

    idx = model.get_idx(model);
    blockIdx = get_block_indices(model, blockNames);
    nParams = numel(blockIdx);
    % Read the blocks to be optimized (active in the energy) from 'additive','g', 'h', 'i', 'j'; blockNames from function build_alternating_schedule.
    solve_additive = any(strcmpi(blockNames, 'additive'));
    solve_h        = any(strcmpi(blockNames, 'h'));
    solve_g        = any(strcmpi(blockNames, 'g'));
    solve_j        = any(strcmpi(blockNames, 'j'));
    solve_i        = any(strcmpi(blockNames, 'i'));
    solve_k        = any(strcmpi(blockNames, 'k'));
    solve_l        = any(strcmpi(blockNames, 'l'));

    nRows = 0;
    for ii = 1:numel(data.modes)
        mode = data.modes{ii};
        Nm = numel(data.(mode).P11);
        nRows = nRows + Nm;
        if ~strcmpi(mode,'SS')
            nRows = nRows + Nm;
        end
    end

    A = zeros(nRows, nParams);
    y = zeros(nRows, 1);

    row = 0;

    for ii = 1:numel(data.modes)
        mode = data.modes{ii};
        D = data.(mode);
        p = prot.(mode);

        Nm = numel(D.P11);
        Pmax = max(abs(D.P11));

        scale1 = sqrt(D.w/Nm/Pmax);        % For experiment-model fitting in the lsq problems.
        scale2 = 100*sqrt(D.w/Nm/Pmax);    % For P22=0 BC in the lsq problems.

        for k = 1:Nm % For all stress-strain data pairs, each corresponding to a raw of A matrix.
            lam = D.lambda(k);
            F = p.Fbar(lam);

            C  = F.' * F;
            I1 = trace(C);
            I2 = 0.5*(I1^2 - trace(C*C));
            J  = sqrt(det(C));
            I1bar = I1 * J^(-2/3);
            I2bar = I2 * J^(-4/3);

            % Fixed current values (evaluate fixed blocks together with optimized blocks)
            h_val = 0; h_der = 0; g_val = 0; g_der = 0;
            j_val = 0; j_der = 0; i_val = 0; i_der = 0;
            k_val = 0; k_der = 0; l_val = 0; l_der = 0;

            if model.use_I1J
                h_val = fnval(model.h_spline, I1bar);
                h_der = fnval(fnder(model.h_spline,1), I1bar);
                g_val = fnval(model.g_spline, J);
                g_der = fnval(fnder(model.g_spline,1), J);
            end

            if model.use_I2J
                j_val = fnval(model.j_spline, I2bar);
                j_der = fnval(fnder(model.j_spline,1), I2bar);
                i_val = fnval(model.i_spline, J);
                i_der = fnval(fnder(model.i_spline,1), J);
            end

            if model.use_I1I2
                k_val = fnval(model.k_spline, I1bar);
                k_der = fnval(fnder(model.k_spline,1), I1bar);
                l_val = fnval(model.l_spline, I2bar);
                l_der = fnval(fnder(model.l_spline,1), I2bar);
            end

            % ---------------- primary row ----------------
            row = row + 1;
            rowA = zeros(1, nParams);
            col = 0;
            % One row with dP11 wrt optimiz; parameters in each column.
            if solve_additive
                if model.use_I1
                    phiI1 = fnval(fnder(model.WI1_spline_dtheta1,1), I1bar);
                    for jcol = 1:numel(phiI1)
                        col = col + 1;
                        Pj = eq_P_from_F(phiI1(jcol), 0.0, 0.0, F);
                        if strcmpi(mode,'SS')
                            rowA(col) = scale1 * Pj(1,2);
                        else
                            rowA(col) = scale1 * Pj(1,1);
                        end
                    end
                end

                if model.use_I2
                    phiI2 = fnval(fnder(model.WI2_spline_dtheta2,1), I2bar);
                    for jcol = 1:numel(phiI2)
                        col = col + 1;
                        Pj = eq_P_from_F(0.0, phiI2(jcol), 0.0, F);
                        if strcmpi(mode,'SS')
                            rowA(col) = scale1 * Pj(1,2);
                        else
                            rowA(col) = scale1 * Pj(1,1);
                        end
                    end
                end

                if model.use_J
                    phiW = fnval(fnder(model.WJ_spline_dtheta,1), J);
                    for jcol = 1:numel(phiW)
                        col = col + 1;
                        Pj = eq_P_from_F(0.0, 0.0, phiW(jcol), F);
                        if strcmpi(mode,'SS')
                            rowA(col) = scale1 * Pj(1,2);
                        else
                            rowA(col) = scale1 * Pj(1,1);
                        end
                    end
                end
            end

            if solve_h
                phiH  = fnval(model.h_spline_dtheta, I1bar);
                dphiH = fnval(fnder(model.h_spline_dtheta,1), I1bar);
                for jcol = 1:numel(phiH)
                    col = col + 1;
                    aj = dphiH(jcol) * g_val;
                    cj = phiH(jcol)  * g_der;
                    Pj = eq_P_from_F(aj, 0.0, cj, F);
                    if strcmpi(mode,'SS')
                        rowA(col) = scale1 * Pj(1,2);
                    else
                        rowA(col) = scale1 * Pj(1,1);
                    end
                end
            end

            if solve_g
                phiG  = fnval(model.g_spline_dtheta, J);
                dphiG = fnval(fnder(model.g_spline_dtheta,1), J);
                for jcol = 1:numel(phiG)
                    col = col + 1;
                    aj = h_der * phiG(jcol);
                    cj = h_val * dphiG(jcol);
                    Pj = eq_P_from_F(aj, 0.0, cj, F);
                    if strcmpi(mode,'SS')
                        rowA(col) = scale1 * Pj(1,2);
                    else
                        rowA(col) = scale1 * Pj(1,1);
                    end
                end
            end

            if solve_j
                phiJ  = fnval(model.j_spline_dtheta, I2bar);
                dphiJ = fnval(fnder(model.j_spline_dtheta,1), I2bar);
                for jcol = 1:numel(phiJ)
                    col = col + 1;
                    bj = dphiJ(jcol) * i_val;
                    cj = phiJ(jcol)  * i_der;
                    Pj = eq_P_from_F(0.0, bj, cj, F);
                    if strcmpi(mode,'SS')
                        rowA(col) = scale1 * Pj(1,2);
                    else
                        rowA(col) = scale1 * Pj(1,1);
                    end
                end
            end

            if solve_i
                phiI  = fnval(model.i_spline_dtheta, J);
                dphiI = fnval(fnder(model.i_spline_dtheta,1), J);
                for jcol = 1:numel(phiI)
                    col = col + 1;
                    bj = j_der * phiI(jcol);
                    cj = j_val * dphiI(jcol);
                    Pj = eq_P_from_F(0.0, bj, cj, F);
                    if strcmpi(mode,'SS')
                        rowA(col) = scale1 * Pj(1,2);
                    else
                        rowA(col) = scale1 * Pj(1,1);
                    end
                end
            end

            if solve_k
                phiK  = fnval(model.k_spline_dtheta, I1bar);
                dphiK = fnval(fnder(model.k_spline_dtheta,1), I1bar);
                for jcol = 1:numel(phiK)
                    col = col + 1;
                    aj = dphiK(jcol) * l_val;   % ∂/∂I1
                    bj = phiK(jcol) * l_der;
                    Pj = eq_P_from_F(aj, bj, 0.0, F);
                    if strcmpi(mode,'SS')
                        rowA(col) = scale1 * Pj(1,2);
                    else
                        rowA(col) = scale1 * Pj(1,1);
                    end
                end
            end

            if solve_l
                phiL  = fnval(model.l_spline_dtheta, I2bar);
                dphiL = fnval(fnder(model.l_spline_dtheta,1), I2bar);
                for jcol = 1:numel(phiL)
                    col = col + 1;
                    aj = phiL(jcol) * k_der;
                    bj = dphiL(jcol) * k_val;   % ∂/∂I2
                    Pj = eq_P_from_F(aj, bj, 0.0, F);
                    if strcmpi(mode,'SS')
                        rowA(col) = scale1 * Pj(1,2);
                    else
                        rowA(col) = scale1 * Pj(1,1);
                    end
                end
            end

            A(row,:) = rowA;
            y(row) = scale1 * D.P11(k); % A*Theta (i.e., P_analytical) should equal experimental D.P11, for optimization parameters Theta.

            % ---------------- P22 row ----------------
            % Additional rows for dP22 wrt optimiz.
            if ~strcmpi(mode,'SS') % P22=0 condition applied only for uniaxial tension and compression.
                row = row + 1;
                rowA22 = zeros(1, nParams);
                col = 0;

                if solve_additive
                    if model.use_I1
                        phiI1 = fnval(fnder(model.WI1_spline_dtheta1,1), I1bar);
                        for jcol = 1:numel(phiI1)
                            col = col + 1;
                            Pj = eq_P_from_F(phiI1(jcol), 0.0, 0.0, F);
                            rowA22(col) = scale2 * Pj(2,2);
                        end
                    end

                    if model.use_I2
                        phiI2 = fnval(fnder(model.WI2_spline_dtheta2,1), I2bar);
                        for jcol = 1:numel(phiI2)
                            col = col + 1;
                            Pj = eq_P_from_F(0.0, phiI2(jcol), 0.0, F);
                            rowA22(col) = scale2 * Pj(2,2);
                        end
                    end

                    if model.use_J
                        phiW = fnval(fnder(model.WJ_spline_dtheta,1), J);
                        for jcol = 1:numel(phiW)
                            col = col + 1;
                            Pj = eq_P_from_F(0.0, 0.0, phiW(jcol), F);
                            rowA22(col) = scale2 * Pj(2,2);
                        end
                    end
                end

                if solve_h
                    phiH  = fnval(model.h_spline_dtheta, I1bar);
                    dphiH = fnval(fnder(model.h_spline_dtheta,1), I1bar);
                    for jcol = 1:numel(phiH)
                        col = col + 1;
                        aj = dphiH(jcol) * g_val;
                        cj = phiH(jcol)  * g_der;
                        Pj = eq_P_from_F(aj, 0.0, cj, F);
                        rowA22(col) = scale2 * Pj(2,2);
                    end
                end

                if solve_g
                    phiG  = fnval(model.g_spline_dtheta, J);
                    dphiG = fnval(fnder(model.g_spline_dtheta,1), J);
                    for jcol = 1:numel(phiG)
                        col = col + 1;
                        aj = h_der * phiG(jcol);
                        cj = h_val * dphiG(jcol);
                        Pj = eq_P_from_F(aj, 0.0, cj, F);
                        rowA22(col) = scale2 * Pj(2,2);
                    end
                end

                if solve_j
                    phiJ  = fnval(model.j_spline_dtheta, I2bar);
                    dphiJ = fnval(fnder(model.j_spline_dtheta,1), I2bar);
                    for jcol = 1:numel(phiJ)
                        col = col + 1;
                        bj = dphiJ(jcol) * i_val;
                        cj = phiJ(jcol)  * i_der;
                        Pj = eq_P_from_F(0.0, bj, cj, F);
                        rowA22(col) = scale2 * Pj(2,2);
                    end
                end

                if solve_i
                    phiI  = fnval(model.i_spline_dtheta, J);
                    dphiI = fnval(fnder(model.i_spline_dtheta,1), J);
                    for jcol = 1:numel(phiI)
                        col = col + 1;
                        bj = j_der * phiI(jcol);
                        cj = j_val * dphiI(jcol);
                        Pj = eq_P_from_F(0.0, bj, cj, F);
                        rowA22(col) = scale2 * Pj(2,2);
                    end
                end

                if solve_k
                    phiK  = fnval(model.k_spline_dtheta, I1bar);
                    dphiK = fnval(fnder(model.k_spline_dtheta,1), I1bar);
                    for jcol = 1:numel(phiK)
                        col = col + 1;
                        aj = dphiK(jcol) * l_val;   % ∂/∂I1
                        bj = phiK(jcol) * l_der;
                        Pj = eq_P_from_F(aj, bj, 0.0, F);
                        rowA22(col) = scale2 * Pj(2,2);
                    end
                end

                if solve_l
                    phiL  = fnval(model.l_spline_dtheta, I2bar);
                    dphiL = fnval(fnder(model.l_spline_dtheta,1), I2bar);
                    for jcol = 1:numel(phiL)
                        col = col + 1;
                        aj = phiL(jcol) * k_der;
                        bj = dphiL(jcol) * k_val;   % ∂/∂I2
                        Pj = eq_P_from_F(aj, bj, 0.0, F);
                        rowA22(col) = scale2 * Pj(2,2);
                    end
                end

                A(row,:) = rowA22;
                y(row) = 0.0; % since P22 should be 0.
            end
            % More rows (as many as experimental stress-strain data pairs) are appended.
        end
    end
end