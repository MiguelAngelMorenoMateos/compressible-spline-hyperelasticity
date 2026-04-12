function [Aineq, bineq] = constraints_hyper(model)
% Linear inequality constraints for:
%   - univariate model: enforce W''(I1) >= 0 and (optionally) W''(I2) >= 0
%   - surface_uv model: enforce What_uu >= 0 and What_vv >= 0 (directional convexity)
%
% Note: full 2D convexity (PSD Hessian) is nonlinear; we cannot enforce it here.

    switch lower(model.type)
        case 'univariate'
            [Aineq, bineq] = constraints_univariate(model);

        case 'surface_uv'
            [Aineq, bineq] = constraints_surface_uv(model);

        otherwise
            Aineq = [];
            bineq = [];
    end

    % Optional: normalize each row (constraints are defined only up to scaling)
    if ~isempty(Aineq)
        s = vecnorm(Aineq, 2, 2);
        s(s==0) = 1;
        Aineq = Aineq ./ s;
    end
end

% =========================== UNIVARIATE ===========================
function [Aineq, bineq] = constraints_univariate(model)
    Aineq = [];
    bineq = [];

    % W(I1) convexity
    if model.use_I1
        B1      = spapi(4, model.I1bar_knots, eye(numel(model.I1bar_knots)));
        B1dd    = fnder(B1, 2);
        A1dd      = -B1dd.coefs';
        b1dd      = zeros(size(A1dd,1), 1);

        % W/(I1) monotonicity: W'(I1) >= 0
        B1d     = fnder(B1, 1);
        A1d     = -B1d.coefs';
        b1d     = zeros(size(A1d,1), 1);    

        A1 = [A1dd; A1d];
        b1 = [ b1dd; b1d ];

        Aineq =  A1;
        bineq = b1;
    end

    if model.use_I2
        B2      = spapi(4, model.I2bar_knots, eye(numel(model.I2bar_knots)));
        B2dd    = fnder(B2, 2);
        A2dd      = -B2dd.coefs';
        b2dd      = zeros(size(A2dd,1), 1);

         % W/(I2) monotonicity: W'(I2) >= 0
        B2d     = fnder(B2, 1);
        A2d     = -B2d.coefs';
        b2d     = zeros(size(A2d,1), 1);

        A2 =[A2dd; A2d];
        b2 = [ b2dd; b2d ]; 

        Aineq = blkdiag(Aineq, A2);
        bineq = [bineq; b2];
    end

    % --- W_vol(J) constraints ---
    if model.use_J
        BJ      = spapi(4, model.J_knots, eye(numel(model.J_knots)));
        BJdd    = fnder(BJ, 2);
        AJdd      = -BJdd.coefs';
        bJdd      = zeros(size(AJdd,1), 1);     

        Aineq = blkdiag(Aineq, AJdd);
        bineq = [bineq; bJdd];
    end


    % --- W(I1,J) = h(I1)*g(J) constraints ---
    if model.use_I1J
        % h(I1) convexity
        Bh = spapi(4, model.h_I1_knots, eye(numel(model.h_I1_knots)));
        Bhdd = fnder(Bh, 2);
        Ahdd = -Bhdd.coefs';
        bhdd = zeros(size(Ahdd,1), 1);

        % h(I1) monotonicity: h'(I1) >= 0
        Bhd     = fnder(Bh, 1);
        Ahd     = -Bhd.coefs';
        bhd     = zeros(size(Ahd,1), 1);

        Ah = [Ahdd; Ahd];
        bh = [ bhdd; bhd ];

        Aineq = blkdiag(Aineq, Ah);
        bineq = [bineq; bh];

        % g(J) convexity
        Bg = spapi(4, model.g_J_knots, eye(numel(model.g_J_knots)));
        Bgdd = fnder(Bg, 2);
        Agdd = -Bgdd.coefs';
        bgdd = zeros(size(Agdd,1), 1);

        % g(J) monotonicity: g'(J) >= 0
        %Bgd     = fnder(Bg, 1);
        %Agd     = -Bgd.coefs';
        %bgd     = zeros(size(Agd,1), 1);

        Ag = Agdd; % enforce only convexity, not monotonicity, for g(J)
        bg = bgdd;

        %Ag = zeros(0, size(Ag,2)); % Uncomment this to disable g(J) constraint
        %bg = zeros(0,1);

        Aineq = blkdiag(Aineq, Ag);
        bineq = [bineq; bg];
    end


    if model.use_I2J
        % --- W(I2,J) = j(I2)*i(J) constraints ---
        % j(I2) convexity
        Bj = spapi(4, model.j_I2_knots, eye(numel(model.j_I2_knots)));
        Bjdd = fnder(Bj, 2);
        Ajdd = -Bjdd.coefs';
        bjdd = zeros(size(Ajdd,1), 1);

        % j(I2) monotonicity: j'(I2) >= 0
        Bjd     = fnder(Bj, 1);
        Ajd     = -Bjd.coefs';
        bjd     = zeros(size(Ajd,1), 1);

        Aj = [Ajdd; Ajd];
        bj = [ bjdd; bjd ];

        Aineq = blkdiag(Aineq, Aj);
        bineq = [bineq; bj];

        % i(J) convexity
        Bi = spapi(4, model.i_J_knots, eye(numel(model.i_J_knots)));
        Bidd = fnder(Bi, 2);
        Aidd = -Bidd.coefs';
        bidd = zeros(size(Aidd,1), 1);

        Ai = Aidd; % enforce only convexity, not monotonicity, for i(J)
        bi = bidd;

        %Ai = zeros(0, size(Ai,2)); % Uncomment this to disable i(J) convexity constraint
        %bi = zeros(0,1);

        Aineq = blkdiag(Aineq, Ai);
        bineq = [bineq; bi];
    end

    if model.use_I1I2
        % --- W(I1,I2) = k(I1)*l(I2) constraints ---
        % k(I1) convexity
        Bk = spapi(4, model.k_I1_knots, eye(numel(model.k_I1_knots)));
        Bkdd = fnder(Bk, 2);
        Akdd = -Bkdd.coefs';
        bkdd = zeros(size(Akdd,1), 1);

        % k(I1) monotonicity: j'(I1) >= 0
        Bkd     = fnder(Bk, 1);
        Akd     = -Bkd.coefs';
        bkd     = zeros(size(Akd,1), 1);

        Ak = [Akdd; Akd];
        bk = [bkdd; bkd ];

        Aineq = blkdiag(Aineq, Ak);
        bineq = [bineq; bk];

        % l(I2) convexity
        Bl = spapi(4, model.l_I2_knots, eye(numel(model.l_I2_knots)));
        Bldd = fnder(Bl, 2);
        Aldd = -Bldd.coefs';
        bldd = zeros(size(Aldd,1), 1);

        % l(I2) monotonicity: j'(I2) >= 0
        Bld     = fnder(Bk, 1);
        Ald     = -Bkd.coefs';
        bld     = zeros(size(Akd,1), 1);

        Al = [Aldd; Ald];
        bl = [bldd; bld ];

        Aineq = blkdiag(Aineq, Al);
        bineq = [bineq; bl];
    end
end


% ============================ SURFACE ============================= 
function [Aineq, bineq] = constraints_surface(model)

    Aineq = [];
    bineq = [];

    % --- W(I1) --- 
    if model.use_I1
        B1      = spapi(4, model.I1bar_knots, eye(numel(model.I1bar_knots)));
        B1dd    = fnder(B1, 2);
        A1dd      = -B1dd.coefs';
        b1dd      = zeros(size(A1dd,1), 1);

        % W/(I1) monotonicity: W'(I1) >= 0
        B1d     = fnder(B1, 1);
        A1d     = -B1d.coefs';
        b1d     = zeros(size(A1d,1), 1);    

        A1 = [A1dd; A1d];
        b1 = [ b1dd; b1d ];

        Aineq =  A1;
        bineq = b1;
    end

    if model.use_I2
        B2      = spapi(4, model.I2bar_knots, eye(numel(model.I2bar_knots)));
        B2dd    = fnder(B2, 2);
        A2dd      = -B2dd.coefs';
        b2dd      = zeros(size(A2dd,1), 1);

        % W/(I2) monotonicity: W'(I2) >= 0
        B2d     = fnder(B2, 1);
        A2d     = -B2d.coefs';
        b2d     = zeros(size(A2d,1), 1);

        A2 =[A2dd; A2d];
        b2 = [ b2dd; b2d ]; 

        Aineq = blkdiag(Aineq, A2);
        bineq = [bineq; b2];
    end

    if model.use_J
        BJ      = spapi(4, model.J_knots, eye(numel(model.J_knots)));
        BJdd    = fnder(BJ, 2);
        AJdd      = -BJdd.coefs';
        bJdd      = zeros(size(AJdd,1), 1);     

        Aineq = blkdiag(Aineq, AJdd);
        bineq = [bineq; bJdd];
    end


    % Surface constraints
    B1d = fnder(model.Wcpl_spline_dtheta, [1 0]); % d^2W/dI1^2 >= 0 (monotonicity in I1 direction)
    B2d = fnder(model.Wcpl_spline_dtheta, [2 0]); % d^2W/dI2^2 >= 0 (convexity in I2 direction)
    BJd = fnder(model.Wcpl_spline_dtheta, [0 1]); % d^2W/dJ^2 >= 0 (monotonicity in J direction)

    % Convert coefficient tensor to linear constraints A*x <= 0
    % coefs is (nParams) x (...) for vector-valued spline
    A1d = -reshape(B1d.coefs, size(B1d.coefs,1), []).';
    A2d = -reshape(B2d.coefs, size(B2d.coefs,1), []).';
    AJd = -reshape(BJd.coefs, size(BJd.coefs,1), []).';

    Acpl = [A1d; A2d; AJd];
    bcpl = zeros(size(Acpl,1), 1);

    Aineq = blkdiag(Aineq, Acpl);
    bineq = [bineq; bcpl];
end
