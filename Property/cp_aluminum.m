function y=cp_aluminum(T)
y=1.199e-5*T.^3-1.516e-2*T.^2+6.512*T-6.285;    %J/kg-K
end
 
% Temperature(K)  Specific heat(J/kg-K)
% 200             790.5
% 250             855.4
% 298.15          897
% 350             930.6
% 400             955.5
% 500             994.8
% 600             1034