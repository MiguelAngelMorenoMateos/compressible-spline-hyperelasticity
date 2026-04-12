function model = eq_params_init_wrapper(I1bar_max, I2bar_max, J_max, J_min, opts)

% Fields returned (common API):
%   model.pack(model) -> x
%   model.unpack(x,model) -> model
%   model.lb() -> lb
%   model.ub() -> ub
%   model.eval_ab(model,F) -> [a,b]
%   model.eval_dabdx(model,F) -> [da_dx, db_dx]
%   model.nParams() -> number of params

    if nargin < 5 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'type') || isempty(opts.type), opts.type = 'multivariate'; end

    switch lower(opts.type)
        case 'multiplicative'
            model = eq_params_init_multiplicative(I1bar_max, I2bar_max, J_max, J_min, opts);
        otherwise
            error('eq_params_init: unknown opts.type = %s', opts.type);
    end

    % Common helper
    model.nParams = @() numel(model.pack(model));
end
