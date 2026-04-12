function model = eq_params_init_wrapper(I1bar_max, I2bar_max, J_max, J_min, opts)
% eq_params_init
% Create equilibrium hyperelastic energy model either as:
%   - separable univariate splines: W(I1)+W(I2)
%   - coupled surface spline: What(u,v) with u=I1, v in [0,1] via admissible mapping
%
% Usage:
%   model = eq_params_init(I1_max, I2_max);  % default = univariate
%   model = eq_params_init(I1_max, I2_max, struct('type','surface'));
%
% Fields returned (common API):
%   model.pack(model) -> x
%   model.unpack(x,model) -> model
%   model.lb() -> lb
%   model.ub() -> ub
%   model.eval_ab(model,F) -> [a,b]
%   model.eval_dabdx(model,F) -> [da_dx, db_dx]
%   model.nParams() -> number of params

    if nargin < 5 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'type') || isempty(opts.type), opts.type = 'univariate'; end % surface or univariate

    switch lower(opts.type)
        case 'univariate'
            model = eq_params_init_univariate(I1bar_max, I2bar_max, J_max, J_min, opts);
        case 'surface'
            model = eq_params_init_surface(I1bar_max, I2bar_max, J_max, J_min, opts);
        case 'multiplicative'
            model = eq_params_init_multiplicative(I1bar_max, I2bar_max, J_max, J_min, opts);
        otherwise
            error('eq_params_init: unknown opts.type = %s', opts.type);
    end

    % Common helper
    model.nParams = @() numel(model.pack(model));
end
