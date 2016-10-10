function [Q_io,h_io]=Annulus(T_i,T_o,r_i,r_o,p_vac)
% W per unit length
kb=1.3806488e-23; %Boltzman constant J/K
T=293; % 20C K
d=4e-10; % diameter of molecule m
% mean free path: kb*T/(sqrt(2)*pi*d^2*p)
l=kb*T/(sqrt(2)*pi*d^2*p_vac*1.01325e+5);
b=1.5833; % b=(2-a)/a*(9*r-5)/2/(r+1) a=1 r=1.4
% thermal conductivity sqrt(kb^3*T/(pi^3*d^4*m))
% p_crit=kb*T/(sqrt(2)*pi*d^2*(r_o-r_i))/1e+5; %bar
% p_crit=3.1*(0.156-0.15)/2/(r_o-r_i)/1e+5; %bar
rho_vac=1.2395; %kg/m3 at 20C
T_avg=(T_i+T_o)/2;
cp=cp_air(T_avg);
mu=mu_air(T_avg);
k=k_air(T_avg);
L_c=2*(log(r_o/r_i))^(4/3)/(r_i^(-3/5)+r_o^(-3/5))^(5/3);
%Prandtl number
Pr_c=cp.*mu./k;  
%Grashof number
Gr=9.8*(max(T_i,T_o)./T_avg-1)*L_c^3./(mu./rho_vac).^2*p_vac^2;
%Rayleigh number
Ra=Gr.*Pr_c;
k_eff=k.*max(1,0.386.*(Pr_c./(0.861+Pr_c)).^0.25.*Ra.^0.25);
h_io=k_eff/(r_i*log(r_o/r_i)+b*l*(r_i/r_o+1))*2*pi*r_i;
Q_io=h_io.*(T_i-T_o);
end
