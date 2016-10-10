function y=mu_air(T)
y=0.1458*T.^1.5./(T+110.5)*1e-5;    %kg/m-s
end

% Temperature(K)  Dynamic viscosity(1e-5 Pa-s or kg/m-s)
% 200             1.329
% 225             1.467
% 250             1.599
% 275             1.725
% 300             1.846
% 325             1.962
% 350             2.075
% 375             2.181
% 400             2.286
% 450             2.485
% 500             2.670