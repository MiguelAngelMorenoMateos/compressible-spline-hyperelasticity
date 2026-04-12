function model = eq_params_init_multiplicative(I1bar_max, I2bar_max, J_max, J_min, opts)

    model.type = 'univariate';
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

    if opts.use_I2J
        model.use_I2J = true;
    else
        model.use_I2J = false;
    end

    if opts.use_I1I2
        model.use_I1I2 = true;
    else
        model.use_I1I2 = false;
    end



    model.k = 4; % cubic splines
    if model.use_I1
        model.n1 = 10; % number of points for I1
    end
    if model.use_I2
        model.n2 = 10;  % number of points for I2
    end
    if model.use_J
        model.nJ = 10; % number of points for J. Important! Set it so that one point of J_knots is exactly at J=1, so that Psi_vol)J=1 = 0 can be enforced.
    end
    if model.use_I1J
        model.nh = 10; % number of points for I1 on spline for h(I1) in multiplicative decomposition. Important!
        model.ng = 10; % number of points for J on spline for g(J) in multiplicative decomposition. Important! Set it so that one point of J_knots is exactly at J=1, so that the constraint g(J=1)=0 can be enforced.
    end
    if model.use_I2J
        model.nj = 10; % number of points for I2 on spline for h(I2) in multiplicative decomposition. Important!
        model.ni = 10; % number of points for J on spline for g(J) in W(I2,J) in multiplicative decomposition. Important! Set it so that one point of J_knots is exactly at J=1, so that the constraint g(J=1)=0 can be enforced.
    end
    if model.use_I1I2
        model.nk = 10;
        model.nl = 10;
    end

    fprintf('Initializing univariate model with use_I1bar_max = %d, use_I2bar_max = %d, J_max = %d, J_min = %d\n', I1bar_max, I2bar_max, J_max, J_min);


    %lam2_free = false;
    %if ~lam2_free
    %    model.I1bar_knots = linspace(3, I1bar_max, model.n1); % sites; NOTE; for lam2 fixed; F known a priori.
    %    model.I2bar_knots = linspace(3, I2bar_max, model.n2); % sites; NOTE; for lam2 fixed; F known a priori.
    %    model.J_knots = linspace(J_min, J_max, model.nJ); % sites; NOTE: for lam2=lam3=1; F known a priori.
    %    %model.J_knots = linspace(0.8, 1, model.nJ); % sites; NOTE: for lam2=lam3=1; F known a priori. For SS
    %else
    %    model.I1bar_knots = linspace(3, 5.5, model.n1); % sites; NOTE: for lam2 free; as variable, F not known a priori.
    %    model.I2bar_knots = linspace(3, 6.5, model.n2); % sites; NOTE: for lam2 free; as variable, F not known a priori.
    %    model.J_knots = linspace(1, 1.6, model.nJ); % manually set the J range to be [0,1.5] since we know the data lives approximately in that range.
    %end

    % --- W(I1) ---
    if model.use_I1
        model.I1bar_knots = linspace(3, I1bar_max, model.n1); % sites; NOTE; for lam2 fixed; F known a priori.
        %model.I1bar_vals  = linspace(0, 0, model.n1); % initial guess
        model.I1bar_vals = 1 * randn(1, model.n1); % Stochastic initialization
        model.WI1_spline         = spapi(model.k, model.I1bar_knots, model.I1bar_vals);
        model.WI1_spline_dtheta1 = spapi(model.k, model.I1bar_knots, eye(model.n1));
    end

    % --- W(I2) ---
    if model.use_I2
        model.I2bar_knots = linspace(3, I2bar_max, model.n2); % sites; NOTE; for lam2 fixed; F known a priori.
        %model.I2bar_vals  = linspace(0, 0, model.n2); % initial guess
        model.I2bar_vals  =  1 * randn(1, model.n2); % Stochastic initialization
        model.WI2_spline         = spapi(model.k, model.I2bar_knots, model.I2bar_vals);
        model.WI2_spline_dtheta2 = spapi(model.k, model.I2bar_knots, eye(model.n2));
    end

    % --- W(J) ---
    if model.use_J
        model.J_knots = linspace(J_min, J_max, model.nJ); % sites; NOTE: for lam2=lam3=1; F known a priori.
        % find index where J = 1 lives in J_knots, and set initial guess to 20 there (roughly matching the scale of WI1,WI2)
        % [~, model.idx_J1] = min(abs(model.J_knots - 1));
        %model.J_vals = linspace(0, 10, model.nJ);
        %model.J_vals = 1 * (model.J_knots - 1) .* (model.J_knots - 1); % initial guess that is 0 at J=0.5 and 10 at J=1.5, roughly matching the scale of WI1,WI2
        model.J_vals = 1 * randn(1, model.nJ); % Stochastic initialization
        model.WJ_spline = spapi(model.k, model.J_knots, model.J_vals);
        model.WJ_spline_dtheta = spapi(model.k, model.J_knots, eye(model.nJ));
    end

    % --- W(I1,J) ---
    if model.use_I1J
        % For multiplicative decomposition, we have W(I1,J) = h(I1)*g(J). We initialize h(I1) and g(J) as follows:
        model.h_I1_knots = linspace(3, I1bar_max, model.nh); % knots for h(I1)
        model.g_J_knots = linspace(J_min, J_max, model.ng); % knots for g(J)
        %model.h_I1_vals = linspace(0, 0, model.nh); % initial guess for h(I1)
        model.h_I1_vals = 1 * randn(1, model.nh); % Stochastic initialization
        %model.g_J_vals = linspace(1, 1, model.ng); % initial guess for g(J); %100 * (model.g_J_knots - 1) .* (model.g_J_knots - 1);
        model.g_J_vals = 1 * randn(1, model.ng); % Stochastic initialization
        model.h_spline = spapi(model.k, model.h_I1_knots, model.h_I1_vals);
        model.g_spline = spapi(model.k, model.g_J_knots, model.g_J_vals);
        model.h_spline_dtheta = spapi(model.k, model.h_I1_knots, eye(model.nh));
        model.g_spline_dtheta = spapi(model.k, model.g_J_knots, eye(model.ng));
    end
    
    % --- W(I2,J) ---
    if model.use_I2J
        model.j_I2_knots = linspace(3, I2bar_max, model.nj); % knots for j(I2)
        model.i_J_knots = linspace(J_min, J_max, model.ni); % knots for i(J) in W(I2,J)
        
        %model.j_I2_vals = linspace(0, 0, model.nj); % initial guess for j(I2)
        model.j_I2_vals = 1 * randn(1, model.nj);
        %model.i_J_vals = linspace(1, 1, model.nj); %1 * (model.i_J_knots - 1) .* (model.i_J_knots - 1); % initial guess for i(J) in W(I2,J)
        model.i_J_vals = 1 * randn(1, model.ni);
        model.j_spline = spapi(model.k, model.j_I2_knots, model.j_I2_vals);
        model.i_spline = spapi(model.k, model.i_J_knots, model.i_J_vals);
        model.j_spline_dtheta = spapi(model.k, model.j_I2_knots, eye(model.nj));
        model.i_spline_dtheta = spapi(model.k, model.i_J_knots, eye(model.ni));
    end

    % --- W(I1,I2) ---
    if model.use_I1I2
        model.k_I1_knots = linspace(3, I1bar_max, model.nk); % knots for j(I2)
        model.l_I2_knots = linspace(3, I2bar_max, model.nl); % knots for i(J) in W(I2,J)
        
        model.k_I1_vals = 1 * randn(1, model.nk);
        %model.k_I1_vals = linspace(0, 0, model.nk);
        model.l_I2_vals = 1 * randn(1, model.nl);
        %model.l_I2_vals = linspace(0, 0, model.nl);
        model.k_spline = spapi(model.k, model.k_I1_knots, model.k_I1_vals);
        model.l_spline = spapi(model.k, model.l_I2_knots, model.l_I2_vals);
        model.k_spline_dtheta = spapi(model.k, model.k_I1_knots, eye(model.nk));
        model.l_spline_dtheta = spapi(model.k, model.l_I2_knots, eye(model.nl));
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
        model.h_spline = make_fnxtr_spline(model.h_spline, model.k);
        model.h_spline_dtheta = make_fnxtr_spline(model.h_spline_dtheta, model.k);
        model.g_spline = make_fnxtr_spline(model.g_spline, model.k);
        model.g_spline_dtheta = make_fnxtr_spline(model.g_spline_dtheta, model.k);
    end
    if model.use_I2J
        model.j_spline = make_fnxtr_spline(model.j_spline, model.k);
        model.i_spline = make_fnxtr_spline(model.i_spline, model.k);
        model.j_spline_dtheta = make_fnxtr_spline(model.j_spline_dtheta, model.k);
        model.i_spline_dtheta = make_fnxtr_spline(model.i_spline_dtheta, model.k);
    end
    if model.use_I1I2
        model.k_spline = make_fnxtr_spline(model.k_spline, model.k);
        model.l_spline = make_fnxtr_spline(model.l_spline, model.k);
        model.k_spline_dtheta = make_fnxtr_spline(model.k_spline_dtheta, model.k);
        model.l_spline_dtheta = make_fnxtr_spline(model.l_spline_dtheta, model.k);
    end

    % --- Packing / Unpacking ---
    model.pack   = @(m) pack_model(m);
    model.unpack = @(x,m) unpack_model(x,m);

    % --- Bounds ---
    model.lb = @() build_lb(model);
    model.ub = @() build_ub(model);

    % --- Indices ---
    model.get_idx = @(m) get_param_indices_impl(m);

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
    %x = [x; model.J_vals(:)];
    if model.use_I1J
        x = [x; model.h_I1_vals(:)];
        x = [x; model.g_J_vals(:)];
    end
    if model.use_I2J
        x = [x; model.j_I2_vals(:)];
        x = [x; model.i_J_vals(:)];
    end
    if model.use_I1I2
        x = [x; model.k_I1_vals(:)];
        x = [x; model.l_I2_vals(:)];
    end
end

function model = unpack_model(x, model)
    n1 = 0;
    n2 = 0;
    if model.use_I1
        n1 = numel(model.I1bar_vals);
        model.I1bar_vals = reshape(x(1:n1), size(model.I1bar_vals));
    else
        n1 = 0;
    end
    if model.use_I2
        n2 = numel(model.I2bar_vals);
        model.I2bar_vals = reshape(x(model.use_I1*n1+1:model.use_I1*n1+n2), size(model.I2bar_vals));
    else
        n2 = 0;
    end
    if model.use_J
        nJ = numel(model.J_vals);
        model.J_vals = reshape(x(model.use_I1*n1+model.use_I2*n2+1:model.use_I1*n1+model.use_I2*n2+nJ), size(model.J_vals));
    else
        nJ = 0;
    end
    if model.use_I1J
        nh = numel(model.h_I1_vals);
        ng = numel(model.g_J_vals);
        model.h_I1_vals = reshape(x(model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+1:model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+nh), size(model.h_I1_vals));
        model.g_J_vals = reshape(x(model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+nh+1:model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+nh+ng), size(model.g_J_vals));
    else
        nh = 0; ng = 0;
    end
    if model.use_I2J
        nj = numel(model.j_I2_vals);
        ni = numel(model.i_J_vals);
        model.j_I2_vals = reshape(x(model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+model.use_I1J*(nh+ng)+1:model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+model.use_I1J*(nh+ng)+nj), size(model.j_I2_vals));
        model.i_J_vals = reshape(x(model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+model.use_I1J*(nh+ng)+nj+1:model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+model.use_I1J*(nh+ng)+nj+ni), size(model.i_J_vals));
    else
        nj = 0; ni = 0;
    end
    if model.use_I1I2
        nk = numel(model.k_I1_vals);
        nl = numel(model.l_I2_vals);
        model.k_I1_vals = reshape(x(model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+model.use_I1J*(nh+ng)+model.use_I2J*(nj+ni)+1:model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+model.use_I1J*(nh+ng)+model.use_I2J*(nj+ni)+nk), size(model.k_I1_vals));
        model.l_I2_vals = reshape(x(model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+model.use_I1J*(nh+ng)+model.use_I2J*(nj+ni)+nk+1:model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+model.use_I1J*(nh+ng)+model.use_I2J*(nj+ni)+nk+nl), size(model.l_I2_vals));
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
        model.h_spline = spapi(model.k, model.h_I1_knots, model.h_I1_vals);
        model.g_spline = spapi(model.k, model.g_J_knots, model.g_J_vals);
    end
    if model.use_I2J
        model.j_spline = spapi(model.k, model.j_I2_knots, model.j_I2_vals);
        model.i_spline = spapi(model.k, model.i_J_knots, model.i_J_vals);
    end
    if model.use_I1I2
        model.k_spline = spapi(model.k, model.k_I1_knots, model.k_I1_vals);
        model.l_spline = spapi(model.k, model.l_I2_knots, model.l_I2_vals);
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
        model.h_spline = make_fnxtr_spline(model.h_spline, model.k);
        model.g_spline = make_fnxtr_spline(model.g_spline, model.k);
    end
    if model.use_I2J
        model.j_spline = make_fnxtr_spline(model.j_spline, model.k);
        model.i_spline = make_fnxtr_spline(model.i_spline, model.k);
    end
    if model.use_I1I2
        model.k_spline = make_fnxtr_spline(model.k_spline, model.k);
        model.l_spline = make_fnxtr_spline(model.l_spline, model.k);
    end
end

function spline = make_fnxtr_spline(spline, order)
    spline = fnxtr(spline, order);
end

function lb = build_lb(model)
    lb_val = Inf;
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
        lb = [lb; 0.0; zeros(numel(model.h_I1_vals)-1,1)];
        lb_g = zeros(numel(model.g_J_vals),1);
        %lb_g = lb_val*ones(numel(model.g_J_vals),1);
        idx_g1 = find(abs(model.g_J_knots - 1) < 1e-6, 1); % Make sure that g_J_knots contains a point at J=1!!! Otherwise no one will be found.
        if isempty(idx_g1)
            error('No knot found at J=1 in g_J_knots. The constraint g(J=1)=0 will not be enforced. Consider adjusting model.g_J_knots to include a point at J=1');
        else
            lb_g(idx_g1) = 1.0;
        end
        lb = [lb; lb_g];
        %lb_g = zeros(numel(model.g_J_vals),1);
        %idx_g1 = find(abs(model.g_J_knots - 1) < 1e-6, 1);
        %lb_g(idx_g1) = 1.0; % g(J=1) = 1.0
        %lb = [lb; lb_g];
    end

    if model.use_I2J
        lb = [lb; 0.0; zeros(numel(model.j_I2_vals)-1,1)];

        lb_i = zeros(numel(model.i_J_vals),1);
        %lb_i = lb_val* ones(numel(model.i_J_vals),1);
        idx_i1 = find(abs(model.i_J_knots - 1) < 1e-6, 1); % Make sure that i_J_knots contains a point at J=1!!! Otherwise no one will be found.
        if isempty(idx_i1)
            error('No knot found at J=1 in i_J_knots. The constraint g(J=1)=0 will not be enforced. Consider adjusting model.i_J_knots to include a point at J=1');
        else
            lb_i(idx_i1) = 1.0;
        end
        lb = [lb; lb_i];
    end

    if model.use_I1I2
        lb = [lb; 0.0; zeros(numel(model.k_I1_vals)-1,1)];
        lb = [lb; 0.0; zeros(numel(model.l_I2_vals)-1,1)];
    end
end

function ub = build_ub(model)
    ub_val = Inf; % Upper bound for all spline values. To help solver, a value smaller than Inf is recommended should not the whole domain of the splines be sampled.
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
        ub = [ub; 0.0; ub_val*ones(numel(model.h_I1_vals)-1,1)];

        ub_g = ub_val*ones(numel(model.g_J_vals),1);
        idx_g1 = find(abs(model.g_J_knots - 1) < 1e-6, 1); % Make sure that g_J_knots contains a point at J=1!!! Otherwise no one will be found.
        if isempty(idx_g1)
            error('No knot found at J=1 in g_J_knots. The constraint g(J=1)=0 will not be enforced. Consider adjusting model.g_J_knots to include a point at J=1');
        else
            ub_g(idx_g1) = 1.0;
        end
        ub = [ub; ub_g];
    end

    if model.use_I2J
        ub = [ub; 0.0; ub_val*ones(numel(model.j_I2_vals)-1,1)];
        ub_i = ub_val*ones(numel(model.i_J_vals),1);
        idx_i1 = find(abs(model.i_J_knots - 1) < 1e-6, 1); % Make sure that i_J_knots contains a point at J=1!!! Otherwise no one will be found.
        if isempty(idx_i1)
            error('No knot found at J=1 in i_J_knots. The constraint g(J=1)=0 will not be enforced. Consider adjusting model.i_J_knots to include a point at J=1');
        else
            ub_i(idx_i1) = 1.0;
        end
        ub = [ub; ub_i];
    end

    if model.use_I1I2
        ub = [ub; 0.0; ub_val*ones(numel(model.k_I1_vals)-1,1)];
        ub = [ub; 0.0; ub_val*ones(numel(model.l_I2_vals)-1,1)];
    end
end

function [a,b,c] = eval_ab_univariate(model, F)
    C  = F.' * F;
    I1 = trace(C);
    I2 = 0.5*(I1^2 - trace(C*C));
    J = sqrt(det(C));
    I1bar = I1 * J^(-2/3);
    I2bar = I2 * J^(-4/3);
    if model.use_I1
        a = fnval(fnder(model.WI1_spline, 1), I1bar); % dPsi(I1)/dI1
    else
        a = 0.0;
    end
    if model.use_I2
        b = fnval(fnder(model.WI2_spline, 1), I2bar); % dPsi(I2)/dI2
    else
        b = 0.0;
    end
    if model.use_J
        c = fnval(fnder(model.WJ_spline, 1), J); % dPsi(J)/dJ
    else
        c = 0.0;
    end

    if model.use_I1J
        a = a + fnval(fnder(model.h_spline, 1), I1bar) * fnval(model.g_spline, J); % dPsi(I1)/dI1 + dPsi(I1,J)/dI1 for multiplicative decomposition 
        c = c + fnval(model.h_spline, I1bar) * fnval(fnder(model.g_spline, 1), J); % dPsi(J)/dJ + dPsi(I1,J)/dJ for multiplicative decomposition
    end
    if model.use_I2J
        b = b + fnval(fnder(model.j_spline, 1), I2bar) * fnval(model.i_spline, J); % dPsi(I2)/dI2 + dPsi(I2,J)/dI2 for multiplicative decomposition
        c = c + fnval(model.j_spline, I2bar) * fnval(fnder(model.i_spline, 1), J); % dPsi(J)/dJ + dPsi(I2,J)/dJ for multiplicative decomposition
    end
    if model.use_I1I2
        a = a + fnval(fnder(model.k_spline, 1), I1bar) * fnval(model.l_spline, I2bar);
        b = b + fnval(model.k_spline, I1bar) * fnval(fnder(model.l_spline, 1), I2bar);
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
        nh = numel(model.h_I1_vals);
        ng = numel(model.g_J_vals);
    else
        nh = 0; ng = 0;
    end
    if model.use_I2J
        nj = numel(model.j_I2_vals);
        ni = numel(model.i_J_vals);
    else
        nj = 0; ni = 0;
    end
    if model.use_I1I2
        nk = numel(model.k_I1_vals);
        nl = numel(model.l_I2_vals);
    else
        nk = 0; nl = 0;
    end
    
    nParams = model.use_I1*n1 + model.use_I2*n2 + model.use_J*nJ;
    if model.use_I1J
        nParams = nParams + nh + ng;
    end
    if model.use_I2J
        nParams = nParams + nj + ni;
    end
    if model.use_I1I2
        nParams = nParams + nk + nl;
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
        % add terms for multiplicative decomposition
        da_dx(model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+1:model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+nh) = fnval(model.g_spline, J) * fnval(fnder(model.h_spline_dtheta, 1), I1bar); 
        da_dx(model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+nh+1:model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+nh+ng) = fnval(fnder(model.h_spline, 1), I1bar) * fnval(model.g_spline_dtheta, J); % product rule

        dc_dx(model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+1:model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+nh) = fnval(fnder(model.g_spline, 1), J) * fnval(model.h_spline_dtheta, I1bar);
        dc_dx(model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+nh+1:model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+nh+ng) = fnval(model.h_spline, I1bar) * fnval(fnder(model.g_spline_dtheta, 1), J); % product rule
    end
    if model.use_I2J
        db_dx(model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+model.use_I1J*(nh+ng)+1:model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+model.use_I1J*(nh+ng)+nj) = fnval(model.i_spline, J) * fnval(fnder(model.j_spline_dtheta, 1), I2bar);
        db_dx(model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+model.use_I1J*(nh+ng)+nj+1:model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+model.use_I1J*(nh+ng)+nj+ni) = fnval(fnder(model.j_spline, 1), I2bar) * fnval(model.i_spline_dtheta, J); % product rule

        dc_dx(model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+model.use_I1J*(nh+ng)+1:model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+model.use_I1J*(nh+ng)+nj) = fnval(fnder(model.i_spline, 1), J) * fnval(model.j_spline_dtheta, I2bar);
        dc_dx(model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+model.use_I1J*(nh+ng)+nj+1:model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+model.use_I1J*(nh+ng)+nj+ni) = fnval(model.j_spline, I2bar) * fnval(fnder(model.i_spline_dtheta, 1), J); % product rule        
    end
    if model.use_I1I2
        da_dx(model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+model.use_I1J*(nh+ng)+model.use_I2J*(nj+ni)+1:model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+model.use_I1J*(nh+ng)+model.use_I2J*(nj+ni)+nk) = fnval(model.l_spline, I2bar) * fnval(fnder(model.k_spline_dtheta, 1), I1bar);
        da_dx(model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+model.use_I1J*(nh+ng)+model.use_I2J*(nj+ni)+nk+1:model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+model.use_I1J*(nh+ng)+model.use_I2J*(nj+ni)+nk+nl) = fnval(fnder(model.k_spline, 1), I1bar) * fnval(model.l_spline_dtheta, I2bar);

        db_dx(model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+model.use_I1J*(nh+ng)+model.use_I2J*(nj+ni)+1:model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+model.use_I1J*(nh+ng)+model.use_I2J*(nj+ni)+nk) = fnval(fnder(model.l_spline,1), I2bar) * fnval(model.k_spline_dtheta, I1bar);
        db_dx(model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+model.use_I1J*(nh+ng)+model.use_I2J*(nj+ni)+nk+1:model.use_I1*n1+model.use_I2*n2+model.use_J*nJ+model.use_I1J*(nh+ng)+model.use_I2J*(nj+ni)+nk+nl) = fnval(fnder(model.l_spline_dtheta,1), I2bar) * fnval(model.k_spline, I1bar);
    end

end


function idx = get_param_indices_impl(model)
    idx = struct();
    pos = 0;

    if model.use_I1
        idx.I1 = pos + (1:model.n1);
        pos = pos + model.n1;
    else
        idx.I1 = [];
    end

    if model.use_I2
        idx.I2 = pos + (1:model.n2);
        pos = pos + model.n2;
    else
        idx.I2 = [];
    end

    if model.use_J
        idx.J = pos + (1:model.nJ);
        pos = pos + model.nJ;
    else
        idx.J = [];
    end

    if model.use_I1J
        idx.h = pos + (1:model.nh);
        pos = pos + model.nh;

        idx.g = pos + (1:model.ng);
        pos = pos + model.ng;
    else
        idx.h = [];
        idx.g = [];
    end

    if model.use_I2J
        idx.j = pos + (1:model.nj);
        pos = pos + model.nj;

        idx.i = pos + (1:model.ni);
        pos = pos + model.ni;
    else
        idx.j = [];
        idx.i = [];
    end

    if model.use_I1I2
        idx.k = pos + (1:model.nk);
        pos = pos + model.nk;

        idx.l = pos + (1:model.nl);
        pos = pos + model.nl;
    else
        idx.k = [];
        idx.l = [];
    end

    idx.additive = [idx.I1, idx.I2, idx.J];
    idx.all = 1:pos;
end