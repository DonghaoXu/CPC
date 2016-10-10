function y=k_aluminum(T)
% y=-1.592e-7*T.^2+1.289e-4*T+0.2131;   %kW/m-K
y=-1.592e-4*T.^2+1.289e-1*T+213.1;  %W/m-K
end

% Temperature(K)  Thermal Conductivity(kW/m-K)
% 250             0.235
% 273             0.236
% 300             0.237
% 350             0.240
% 400             0.240
% 500             0.237
% 600             0.232
% 700             0.226