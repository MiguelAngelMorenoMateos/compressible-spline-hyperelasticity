function D = load_foam_turbo(mode)
% mode = 'UT' | 'UC' | 'SS' 
mode = upper(mode);

fname = fullfile('data','experimental','foam_turbo', sprintf('foam_%s.csv', mode));
T = readtable(fname);

D = struct();
D.lambda = T.lambda;
D.P11    = T.P_MPa;
D.dt     = [];   % no time data for Treloar
D.w      = 1.0;  % optional weight
end
