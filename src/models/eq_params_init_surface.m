function model = eq_params_init_surface(I1bar_max, I2bar_max, J_max, J_min, opts)

    model.type = 'surface';
    if ~isfield(opts,'use_I2') || isempty(opts.use_I2)
        model.use_I2 = true;
    else
        model.use_I2 = logical(opts.use_I2);
    end

    if opts.use_I1
        model.use_I1 = true;
    else
        model.use_I1 = false;
    end

    if opts.use_J
        model.use_J = true;
    else
        model.use_J = false;
    end

    if opts.use_I1J
        model.use_I1J = true;
    else
        model.use_I1J = false;
    end



    model.k = 4; % cubic splines
    model.n1 = 3; % number of points for I1
    model.n2 = 5;  % number of points for I2
    model.nJ = 10; % number of points for J. Important! Set it so that one point of J_knots is exactly at J=1, so that Psi_vol)J=1 = 0 can be enforced.
    model.nh = 10; % number of points for I1 on surface spline for W(I1,J). Important!
    model.ng = 10; % number of points for J on surface spline for W(I1,J). Important!

    fprintf('Initializing surface model with use_I1bar_max = %d, use_I2bar_max = %d, J_max = %d, J_min = %d\n', I1bar_max, I2bar_max, J_max, J_min);


    lam2_free = false;
    if ~lam2_free
        model.I1bar_knots = linspace(3, I1bar_max, model.n1); % sites; NOTE; for lam2 fixed; F known a priori.
        model.I2bar_knots = linspace(3, I2bar_max, model.n2); % sites; NOTE; for lam2 fixed; F known a priori.
        model.J_knots = linspace(J_min, J_max, model.nJ); % sites; NOTE: for lam2=lam3=1; F known a priori.
        %model.J_knots = linspace(0.8, 1, model.nJ); % sites; NOTE: for lam2=lam3=1; F known a priori. For SS
    else
        model.I1bar_knots = linspace(3, 5.5, model.n1); % sites; NOTE: for lam2 free; as variable, F not known a priori.
        model.I2bar_knots = linspace(3, 6.5, model.n2); % sites; NOTE: for lam2 free; as variable, F not known a priori.
        model.J_knots = linspace(1, 1.6, model.nJ); % manually set the J range to be [0,1.5] since we know the data lives approximately in that range.
    end

    % --- W(I1) ---
    if model.use_I1
        model.I1bar_vals  = linspace(0, 0, model.n1); % initial guess
        model.WI1_spline         = spapi(model.k, model.I1bar_knots, model.I1bar_vals);
        model.WI1_spline_dtheta1 = spapi(model.k, model.I1bar_knots, eye(model.n1));
    end

    % --- W(I2) ---
    if model.use_I2
        model.I2bar_vals  = linspace(0, 20, model.n2); % initial guess
        model.WI2_spline         = spapi(model.k, model.I2bar_knots, model.I2bar_vals);
        model.WI2_spline_dtheta2 = spapi(model.k, model.I2bar_knots, eye(model.n2));
    end

    % --- W(J) ---
    if model.use_J
        model.J_vals = 10 * (model.J_knots - 1) .* (model.J_knots - 1); % initial guess that is 0 at J=0.5 and 10 at J=1.5, roughly matching the scale of WI1,WI2
        model.WJ_spline = spapi(model.k, model.J_knots, model.J_vals);
        model.WJ_spline_dtheta = spapi(model.k, model.J_knots, eye(model.nJ));
    end

    % --- W(I1,J) ---
    if model.use_I1J
        nh = model.nh; ng = model.ng;
        model.h_sites     = linspace(3, I1bar_max, nh);
        model.g_sites     = linspace(J_min, J_max, ng);

        % initial guess
        W1              = spapi(model.k, model.h_sites, linspace(0, 20, nh));
        W2              = spapi(model.k, model.g_sites, linspace(0, 40, ng));
        model.Wcpl_vals = zeros(nh, ng);
        [uu, vv]    = ndgrid(model.h_sites, model.g_sites);
        for j = 1:ng
            u = uu(:,j);
            v = vv(:,j);
            model.Wcpl_vals(:,j) = fnval(W1, u) + fnval(W2, v);
        end

        % surface spline
        model.Wcpl_spline = spapi({model.k, model.k} , {model.h_sites, model.g_sites}, model.Wcpl_vals);

        % plot surface spline to check
        %fnplt(model.Wcpl_spline); 
        %abcd
    
        % sensitivity surface wrt parameters x = Wcpl_vals(:)
        Eval    = reshape(eye(nh*ng), [nh*ng, nh, ng]);
        model.Wcpl_spline_dtheta = spapi({model.k, model.k}, {model.h_sites, model.g_sites}, Eval);
    end


    % --- Make fnxtr splines for extrapolation ---
    if model.use_I1
        model.WI1_spline = make_fnxtr_spline(model.WI1_spline, model.k);
        model.WI1_spline_dtheta1 = make_fnxtr_spline(model.WI1_spline_dtheta1, model.k);
    end
    if model.use_I2
        model.WI2_spline = make_fnxtr_spline(model.WI2_spline, model.k);
        model.WI2_spline_dtheta2 = make_fnxtr_spline(model.WI2_spline_dtheta2, model.k);
    end
    if model.use_J
        model.WJ_spline = make_fnxtr_spline(model.WJ_spline, model.k);
        model.WJ_spline_dtheta = make_fnxtr_spline(model.WJ_spline_dtheta, model.k);
    end
    if model.use_I1J
        model.Wcpl_spline = make_fnxtr_spline(model.Wcpl_spline, model.k);
        model.Wcpl_spline_dtheta = make_fnxtr_spline(model.Wcpl_spline_dtheta, model.k);
    end
    
    
    % --- Packing / Unpacking ---
    model.pack   = @(m) pack_model(m);
    model.unpack = @(x,m) unpack_model(x,m);

    % --- Bounds ---
    model.lb = @() build_lb(model);
    model.ub = @() build_ub(model);

    % --- Formulation-agnostic wrappers ---
    model.eval_ab    = @(m,F) eval_ab_univariate(m,F);
    model.eval_dabdx = @(m,F) eval_dabdx_univariate(m,F);
end

function x = pack_model(model)
    x = [];
    if model.use_I1
        x = model.I1bar_vals(:);
    end
    if model.use_I2
        x = [x; model.I2bar_vals(:)];
    end
    if model.use_J
        x = [x; model.J_vals(:)];
    end
    if model.use_I1J
        x = [x; model.Wcpl_vals(:)];
    end
end

function model = unpack_model(x, model)
    n1 = 0;
    n2 = 0;
    nJ = 0;
    if model.use_I1
        n1 = numel(model.I1bar_vals);
        model.I1bar_vals = reshape(x(1:n1), size(model.I1bar_vals));
    end
    if model.use_I2
        n2 = numel(model.I2bar_vals);
        model.I2bar_vals = reshape(x(model.use_I1*n1+1:model.use_I1*n1+n2), size(model.I2bar_vals));
    end
    if model.use_J
        nJ = numel(model.J_vals);
        model.J_vals = reshape(x(model.use_I1*n1+model.use_I2*n2+1:model.use_I1*n1+model.use_I2*n2+nJ), size(model.J_vals));
    end
    if model.use_I1J
        nh = model.nh;
        ng = model.ng;
        model.Wcpl_vals = reshape(x(model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+1:end), [nh, ng]);
    end

    if model.use_I1
        model.WI1_spline = spapi(model.k, model.I1bar_knots, model.I1bar_vals);
    end
    if model.use_I2
        model.WI2_spline = spapi(model.k, model.I2bar_knots, model.I2bar_vals);
    end
    if model.use_J
        model.WJ_spline = spapi(model.k, model.J_knots, model.J_vals);
    end

    if model.use_I1J
        model.Wcpl_spline = spapi({model.k, model.k} , {model.h_sites, model.g_sites}, model.Wcpl_vals);
    end

    % --- Make fnxtr splines for extrapolation ---
    if model.use_I1
        model.WI1_spline = make_fnxtr_spline(model.WI1_spline, model.k);
    end
    if model.use_I2
        model.WI2_spline = make_fnxtr_spline(model.WI2_spline, model.k);
    end
    if model.use_J
        model.WJ_spline = make_fnxtr_spline(model.WJ_spline, model.k);
    end
    if model.use_I1J
        model.Wcpl_spline = make_fnxtr_spline(model.Wcpl_spline, model.k);
    end
end



function spline = make_fnxtr_spline(spline, order)
    spline = fnxtr(spline, order);
end



function lb = build_lb(model)
    lb_val = -Inf;
    lb = [];
    if model.use_I1
        lb = [lb; 0.0; zeros(numel(model.I1bar_vals)-1,1)];
    end
    if model.use_I2
        lb = [lb; 0.0; zeros(numel(model.I2bar_vals)-1,1)];
    end
    if model.use_J
        lb = [lb; 0.0; zeros(numel(model.J_vals)-1,1)];
    end

    if model.use_I1J
        lb = [lb; zeros(numel(model.Wcpl_vals(:)),1)];
    end
end

function ub = build_ub(model)
    ub_val = 5.0; % Upper bound for all spline values. To help solver, a value smaller than Inf is recommended should not the whole domain of the splines be sampled.
    ub = [];
    if model.use_I1
        ub = [ub; 0.0; ub_val*ones(numel(model.I1bar_vals)-1,1)];
    end
    if model.use_I2
        ub = [ub; 0.0; ub_val*ones(numel(model.I2bar_vals)-1,1)];
    end
    if model.use_J
        ub_J = ub_val*ones(numel(model.J_vals),1);
        % with idx_J1 the position where J_knots = 1.
        idx_J1 = find(abs(model.J_knots - 1) < 1e-6, 1); % Make sure that J_knots contains a point at J=1!!! Otherwise no one will be found.
        % If no value for idx_J1 is found, write a worning on the screen to say that the constraint Psi_vol(J=1)=0 will not be enforced, and stop all running codes
        if isempty(idx_J1)
            error('No knot found at J=1. The constraint Psi_vol(J=1)=0 will not be enforced. Consider adjusting model.J_knots to include a point at J=1');
        else
            ub_J(idx_J1) = 0.0;
        end
        ub = [ub; ub_J];
    end
    
    if model.use_I1J
        ub = [ub; ub_val*ones(numel(model.Wcpl_vals(:)),1)];
    end
end

function [a,b,c] = eval_ab_univariate(model, F)
    C  = F.' * F;
    I1 = trace(C);
    I2 = 0.5*(I1^2 - trace(C*C));
    J = sqrt(det(C));
    I1bar = I1 * J^(-2/3);
    I2bar = I2 * J^(-4/3);

    a = 0.0; b = 0.0; c = 0.0;
    if model.use_I1
        a = fnval(fnder(model.WI1_spline, 1), I1bar); % dPsi(I1)/dI1
    end
    if model.use_I2
        b = fnval(fnder(model.WI2_spline, 1), I2bar); % dPsi(I2)/dI2
    end
    if model.use_J
        c = fnval(fnder(model.WJ_spline, 1), J); % dPsi(J)/dJ
    end

    if model.use_I1J
        a = a + fnval(fnder(model.Wcpl_spline, [1 0]), [I1bar; J]); % dPsi(I1,J)/dI1
        c = c + fnval(fnder(model.Wcpl_spline, [0 1]), [I1bar; J]); % dPsi(I1,J)/dJ
    end
end



function [da_dx, db_dx, dc_dx] = eval_dabdx_univariate(model, F)  % Compute second derivatives for Jacobian in optimization. 
    C  = F.' * F;
    I1 = trace(C);
    I2 = 0.5*(I1^2 - trace(C*C));
    J = sqrt(det(C));
    I1bar = I1 * J^(-2/3);
    I2bar = I2 * J^(-4/3);

    if model.use_I1
        n1 = model.n1;
        da_dx1 = fnval(fnder(model.WI1_spline_dtheta1, 1), I1bar);
    else
        n1 = 0;
        db_dx1 = zeros(0,1);
    end
    if model.use_I2
        n2 = model.n2;
        db_dx2 = fnval(fnder(model.WI2_spline_dtheta2, 1), I2bar);
    else
        n2 = 0;
        db_dx2 = zeros(0,1);
    end
    if model.use_J
        dc_dx3 = fnval(fnder(model.WJ_spline_dtheta, 1), J);
        nJ = model.nJ;
    else
        dc_dx3 = 0.0;
        nJ = 0;
    end

    if model.use_I1J
        nh = model.nh;
        ng = model.ng;
    end
    
    nParams = model.use_I1*n1 + model.use_I2*n2 + model.use_J*nJ;
    if model.use_I1J
        nParams = nParams + nh * ng;
    end

    da_dx = zeros(nParams,1);
    db_dx = zeros(nParams,1);
    dc_dx = zeros(nParams,1);

    if model.use_I1
        da_dx(1:n1) = da_dx1;
    end
    if model.use_I2
        db_dx(model.use_I1*n1+1:model.use_I1*n1+n2) = db_dx2;
    end
    if model.use_J
        dc_dx(model.use_I1*n1+model.use_I2*n2+1:model.use_I1*n1+model.use_I2*n2+nJ) = dc_dx3;
    end
    
    if model.use_I1J
        da_dx(model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+1:model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+nh*ng) = fnval(fnder(model.Wcpl_spline_dtheta, [1 0]), [I1bar; J]);
        dc_dx(model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+1:model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+nh*ng) = fnval(fnder(model.Wcpl_spline_dtheta, [0 1]), [I1bar; J]);
    end
end
