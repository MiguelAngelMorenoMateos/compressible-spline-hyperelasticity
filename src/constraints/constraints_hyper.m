function [Aineq, bineq] = constraints_hyper(model)
% Linear inequality constraints for:
%   - univariate model: W(I1,I2,J) = W1(I1) + W2(I2) + W_vol(J)
%   - separable model: W(I1,I2,J) = W1(I1) + W2(I2) + W_vol(J) + h(I1)*g(J)
%   - general model: W(I1,I2,J) = W1(I1) + W2(I2) + W_vol(J) + h(I1)*g(J) + j(I2)*i(J) + k(I1)*l(I2)

    switch lower(model.type)
        case 'univariate'
            [Aineq, bineq] = constraints_univariate(model);

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

        % Two inequality constraints:
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

        Ag = Agdd; % enforce only convexity, not monotonicity, for g(J)
        bg = bgdd;

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

        % k(I1) monotonicity: k'(I1) >= 0
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

        % l(I2) monotonicity: l'(I2) >= 0
        Bld     = fnder(Bk, 1);
        Ald     = -Bkd.coefs';
        bld     = zeros(size(Akd,1), 1);

        Al = [Aldd; Ald];
        bl = [bldd; bld ];

        Aineq = blkdiag(Aineq, Al);
        bineq = [bineq; bl];
    end
end