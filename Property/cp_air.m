function y=cp_air(T)
% y=3.451e-7*T.^2-1.529e-4*T+1.02;    %kJ/kg-K
y=3.451e-4*T.^2-1.529e-1*T+1020;    %J/kg-K
end

% Temperature(K)  Specific heat(kJ/kg-K)
% 200             1.0025
% 225             1.0027
% 250             1.0031
% 275             1.0038
% 300             1.0049
% 325             1.0063
% 350             1.0082
% 375             1.0106
% 400             1.0135
% 450             1.0206
% 500             1.0295