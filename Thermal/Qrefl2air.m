function [Qi,hi,Qo,ho]=Qrefl2air(Vel,T_refl,T_amb,height,width,rho_air_a,mu_air_a,cp_air_a,k_air_a,lcpc)
% W per unit length
% global rho_air_a mu_air_a k_air_a cp_air_a lcpc
% In the cavity
Gr_r=9.8*max(1/300,(T_refl/T_amb-1))*height^3/(mu_air_a/rho_air_a)^2;
Rei=max(Vel*width*rho_air_a/mu_air_a,1);
Nui=7.96*Rei^0.18*(Gr_r/Rei^2).^(-0.02)*(height/width)^1.1;
hi=Nui*k_air_a/height; %W/m2-K
Qi=lcpc*hi.*(T_refl-T_amb); %W/m

% Backside
Pr=cp_air_a*mu_air_a/k_air_a;
Reo=max(Vel*2*height*rho_air_a/mu_air_a,1);
Nuo=0.153*Reo^0.638*Pr^(1/3);
ho=Nuo*k_air_a/2/height;
ho=ho*ones(size(T_refl));
Qo=lcpc*ho.*(T_refl-T_amb);
end