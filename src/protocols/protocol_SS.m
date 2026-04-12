% src/protocols/protocol_SS.m
function protocol = protocol_SS()
    protocol.name = "simple_shear";

    % returns deformation gradient F for given stretch ratio lambda
    stretch1 = 0.8; % Constant pre-stretch in the axial direction.
    protocol.Fbar = @(lambda) [stretch1, lambda, 0; 0, 1, 0; 0, 0, 1];

    % return the maximal value of the first and second invariants for given lambda_max
    J = stretch1;

    protocol.I1bar_max = @(lambda_max) J^(-2/3)*(stretch1^2+lambda_max^2+2);
    protocol.I2bar_max = @(lambda_max) J^(-4/3)*(2*stretch1^2+lambda_max^2+1);

end
