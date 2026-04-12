function norms = sparsity_norms(model)
    norms = [];
    if model.use_I1
        W_I1_ctrl   = fnval(model.WI1_spline, model.I1bar_knots);
        % Calculate L2 norm of the coefficients of the splines
        L2_WI1_norm = norm(W_I1_ctrl)/model.n1; % normalize by number of coefficients to make it comparable across different spline sizes
        norms = [norms, L2_WI1_norm];
        fprintf('L2 norm of W(I1) coefficients: %.4f\n', L2_WI1_norm);
    end
    if model.use_I2
        W_I2_ctrl   = fnval(model.WI2_spline, model.I2bar_knots);
        L2_WI2_norm = norm(W_I2_ctrl)/model.n2; % normalize by number of coefficients to make it comparable across different spline sizes
        norms = [norms, L2_WI2_norm];
        fprintf('L2 norm of W(I2) coefficients: %.4f\n', L2_WI2_norm);
    end
    if model.use_J
        W_J_ctrl   = fnval(model.WJ_spline, model.J_knots);
        L2_WJ_norm = norm(W_J_ctrl)/model.nJ; % normalize by number of coefficients to make it comparable across different spline sizes
        norms = [norms, L2_WJ_norm];
        fprintf('L2 norm of W(J) coefficients: %.4f\n', L2_WJ_norm);
    end 
    if model.use_I1J
        h_ctrl = fnval(model.h_spline, model.h_I1_knots);
        g_ctrl = fnval(model.g_spline, model.g_J_knots);
        L2_h_norm = norm(h_ctrl)/model.nh; % normalize by number of coefficients to make it comparable across different spline sizes
        L2_g_norm = norm(g_ctrl)/model.ng; % normalize by number of coefficients to make it comparable across different spline sizes
        L2_hg_norm = norm([L2_h_norm, L2_g_norm])/2; % combined norm for the coupling terms
        norms = [norms, L2_hg_norm];
        fprintf('L2 norm of combined h(I1)g(J) coefficients: %.4f\n', L2_hg_norm);
    end
    if model.use_I2J
        j_ctrl = fnval(model.j_spline, model.j_I2_knots);
        i_ctrl = fnval(model.i_spline, model.i_J_knots);
        L2_j_norm = norm(j_ctrl)/model.nj; % normalize by number of coefficients to make it comparable across different spline sizes
        L2_i_norm = norm(i_ctrl)/model.ni; % normalize by number of coefficients to make it comparable across different spline sizes
        L2_ji_norm = norm([L2_j_norm, L2_i_norm])/2; % combined norm for the coupling terms
        norms = [norms, L2_ji_norm];
        fprintf('L2 norm of combined j(I2)i(J) coefficients: %.4f\n', L2_ji_norm);
    end
    if model.use_I1I2
        k_ctrl = fnval(model.k_spline, model.k_I1_knots);
        l_ctrl = fnval(model.l_spline, model.l_I2_knots);
        L2_k_norm = norm(k_ctrl)/model.nk; % normalize by number of coefficients to make it comparable across different spline sizes
        L2_l_norm = norm(l_ctrl)/model.nl; % normalize by number of coefficients to make it comparable across different spline sizes
        L2_kl_norm = norm([L2_k_norm, L2_l_norm])/2; % combined norm for the coupling terms
        norms = [norms, L2_kl_norm];
        fprintf('L2 norm of k(I1)l(I2) coefficients: %.4f\n', L2_kl_norm);
    end
end