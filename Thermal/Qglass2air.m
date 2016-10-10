function [Q,h]=Qglass2air(Vel,T_amb,T_g,r_g,height,width,rho_air_a,mu_air_a,cp_air_a,k_air_a)
% W per unit length
% global rho_air_a mu_air_a cp_air_a k_air_a
Re=Vel*2*r_g*rho_air_a/mu_air_a;
Gr=9.8*(T_g./T_amb-1)*height^3/(mu_air_a/rho_air_a)^2;
Pr=cp_air_a*mu_air_a/k_air_a;
if Re<1
    Nu=0.44*(2*height/width)^(1/6.5)*(Gr*Pr).^0.25;
    h=Nu*k_air_a/height; %W/m2-K
else
    Pr_s=cp_air(T_g).*mu_air(T_g)./k_air(T_g);
    [C,m]=coefCm(Re);
    Nu=C*Re^m*Pr^0.37*(Pr./Pr_s).^0.25;
    h=Nu*k_air_a/r_g/2; %W/m2-K
end
Q=h*2*pi*r_g.*(T_g-T_amb); %W/m
end

function [C,m]=coefCm(Re)
if Re<40
    C=0.75;
    m=0.4;
elseif Re<1000
    C=0.51;
    m=0.5;
elseif Re<200000
    C=0.26;
    m=0.6;
elseif Re<1000000
    C=0.076;
    m=0.7;
end
end