function run_fit_eq()

    % ----------------------------
    % 1) Load or synthesize data
    % ----------------------------
    type = 'foam_leap'; % 'foam_leap' | 'foam_turbo'
    modes = {'UT','UC','SS'}; % 'UT' | 'UC' | 'SS'
    switch type
        case 'foam_leap'
            data = get_data(type, [], modes);
        case 'foam_turbo'
            data = get_data(type, [], modes);
        otherwise
            error('Unknown data type');
    end
    % print data summary
    fprintf('Loaded data type: %s\n', type);
    fprintf('Modes: '); fprintf('%s ', modes{:}); fprintf('\n');
    for mode = modes
        mode = mode{1};
        D = data.(mode);
        fprintf(' Mode: %s, number of data points: %d\n', mode, numel(D.P11));
    end

    % ----------------------------
    % 2) Model + protocol handles
    % ----------------------------
    prot = protocol_map();

    I1bar_max = -inf; I2bar_max = -inf;
    J_max = 1.0; J_min = 0.9; % errors with +- Inf and SS only
    for i = 1:numel(data.modes)
        mode = data.modes{i};
        p = prot.(mode);
        lam_max = max(data.(mode).lambda);
        % For compression, the maximal stretch is the minimal lambda
        if strcmpi(mode,'UC')
            lam_max = min(data.(mode).lambda);
        end
        I1bar_max = max(I1bar_max, p.I1bar_max(lam_max));
        I2bar_max = max(I2bar_max, p.I2bar_max(lam_max));
        % if not:
        if strcmpi(mode, 'UC') || strcmpi(mode, 'UT')
            J_max = max(J_max, max(data.(mode).lambda));
            J_min = min(J_min, min(data.(mode).lambda));
        end
    end

    c = 1.0; % safety factor for initial guess
    opts.type = 'multiplicative'; % 'multiplicative' or 'surface'
    opts.use_I1 = true;
    opts.use_I2 = true;
    opts.use_J = true;
    opts.use_I1J = true;
    opts.use_I2J = false;
    opts.use_I1I2 = false;
    model = eq_params_init_wrapper(I1bar_max * c, I2bar_max * c, J_max * c, J_min * c, opts); % initial guess with slightly enlarged domain

    % ----------------------------
    % 3) Pack parameters for solver (lsqlin)
    % ----------------------------
    x0 = model.pack(model);
    lb = model.lb();
    ub = model.ub();
    fprintf('Number of parameters to fit: %d\n', numel(x0));
    fprintf('  lb: '); fprintf('%.3f ', lb); fprintf('\n');
    fprintf('  ub: '); fprintf('%.3f ', ub); fprintf('\n');  
    fprintf('  x0: '); fprintf('%.3f ', x0); fprintf('\n');


    % Linear inequality constraints: 
    [Aineq, bineq] = constraints_hyper(model); 
    fprintf('Number of linear constraints: %d\n', size(Aineq,1)); 

    % dPsi_vol/dJ (J=1) = 0 as an equality constraint to ensure zero stress at underformed reference state.
    Aeq = []; beq = [];
    pass_Aeq = true;

    if pass_Aeq
        % Two equality constraints:
        % (1) dPsi_vol/dJ (J=1) = 0
        Aeq = zeros(0, numel(x0)); beq = zeros(0,1);
        if model.use_J
            idx_J = (1:model.nJ); % indices of J spline values in x
            if model.use_I1, idx_J=idx_J+model.n1; end
            if model.use_I2, idx_J=idx_J+model.n2; end 
            Aeq(end+1, idx_J) = fnval(fnder(model.WJ_spline_dtheta, 1), 1);
            beq(end+1, 1) = 0;
            fprintf('Added equality constraint for zero stress at zero stretch: dPsi_vol/dJ (J=1) = 0\n');
        end
        
        % normalize each row of Aeq to have unit norm (since constraints are defined only up to scaling)
        s = vecnorm(Aeq, 2, 2); %p = 2 → Euclidean norm; dim = 2 → operate along columns → meaning compute row-wise norms
        Aeq = Aeq ./ s;
        beq = beq ./ s;
    end

    % Solver options
    tol = 1e-10;
    
    opts = optimoptions('lsqlin', ...
        'Algorithm','active-set', ... % 'interior-point' or 'active-set'
        'Display','iter', ...
        'OptimalityTolerance',tol, ...
        'FunctionTolerance',tol, ...
        'StepTolerance',tol);


    % RUN OPTIMIZATION: LINEAR ALTERNATING SCHEME
    schedule = build_alternating_schedule(model);
    maxAltIter = 300;
    prevObj = inf;

    for it = 1:maxAltIter
        fprintf('\n=== Alternating iteration %d ===\n', it);

        for s = 1:numel(schedule)
            blockNames = schedule{s};
            [model, out_blk] = solve_block_lsqlin(model, prot, data, Aineq, bineq, Aeq, beq, x0, opts, blockNames);

            if iscell(blockNames)
                fprintf('Solved blocks: %s\n', strjoin(blockNames, ', '));
            else
                fprintf('Solved block: %s\n', blockNames);
            end
        end

        x_full = model.pack(model);
        r_full = loss_function(x_full, model, prot, data);
        obj = sqrt(sum(r_full.^2));
        fprintf('  full objective = %.6e\n', obj);

        % Print the norms of the spline coefficients to monitor sparsity during optimization.
        norms_ = sparsity_norms(model);
        fprintf(' Sparsity norms of each spline: \n');
        fprintf('  %s\n', num2str(norms_));
        
        if abs(prevObj - obj) / max(1,abs(prevObj)) < 5e-8
            fprintf('Alternating loop converged.\n');
            break;
        end
        prevObj = obj;
    end
    x_opt = model.pack(model);
    residual = loss_function(x_opt, model, prot, data); % Get final residual and Jacobian at the optimal solution.

    fprintf('Fitted parameters:\n');
    fprintf(' x = '); fprintf('%.4f ', x_opt); fprintf('\n');

    % ----------------------------
    % 5) Save results
    % ----------------------------
    model = model.unpack(x_opt, model);
    result = struct();
    result.x_opt = x_opt;
    result.residual = residual;
    result.norms = norms_;
    result.model = model;
    result.prot = prot;
    result.data = data;
    save('fit_result_eq.mat', 'result');
end