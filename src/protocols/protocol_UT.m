% src/protocols/protocol_UT.m
function protocol = protocol_UT()
    protocol.name = "uniax";

    % returns deformation gradient F for given stretch ratio lambda
    protocol.Fbar = @(lambda) diag([lambda, 1 / 1, 1 / 1]);

    % return the maximal value of the first and second invariants for given lambda_max
    protocol.I1bar_max = @(lambda_max) lambda_max^(-2/3) * (lambda_max^2 + 2);
    protocol.I2bar_max = @(lambda_max) lambda_max^(-4/3) * (2 * lambda_max^2 + 1);
    protocol.J_max = @(lambda_max) lambda_max * (1 / 1) * (1 / 1); % for lam2=lam3=1; F known a priori.
    protocol.J_min = @(lambda_min) lambda_min * (1 / 1) * (1 / 1); % for lam2=lam3=1; F known a priori.
end
