classdef convection
    methods (Static)
        function [Q,h]=Qglass2air(Vel,T_amb,T_g,r_g,height,width,varargin)
            % W per unit length
            in=[varargin(:);cell(4,1)]; in=in(1:4);
            if isempty(in{1})
                in{1}=air();
                in{1}.setdefault(T_amb);
            end
            if isa(in{1},'air')
                rho_air_a=in{1}.rho_def;
                mu_air_a=in{1}.mu_def;
                cp_air_a=in{1}.cp_def;
                k_air_a=in{1}.k_def;
            else
                rho_air_a=in{1};
                mu_air_a=in{2};
                cp_air_a=in{3};
                k_air_a=in{4};
            end
            Re=Vel*2*r_g*rho_air_a/mu_air_a;
            Gr=9.8*(T_g./T_amb-1)*height^3/(mu_air_a/rho_air_a)^2;
            Pr=cp_air_a*mu_air_a/k_air_a;
            if Re<1
                Nu=0.44*(2*height/width)^(1/6.5)*(Gr*Pr).^0.25;
                h=Nu*k_air_a/height; %W/m2-K
            else
                Pr_s=cp_air(T_g).*mu_air(T_g)./k_air(T_g);
                [C,m]=convection.coefCm(Re);
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
        
        function [Qi,hi,Qo,ho]=Qrefl2air(Vel,T_refl,T_amb,height,width,varargin)
            % W per unit length
            in=[varargin(:);cell(5,1)]; in=in(1:5);
            lcpc=in{1};
            if isempty(in{2})
                in{2}=air();
                in{2}.setdefault(T_amb);
            end    
            if isa(in{2},'air')
                rho_air_a=in{2}.rho_def;
                mu_air_a=in{2}.mu_def;
                cp_air_a=in{2}.cp_def;
                k_air_a=in{2}.k_def;
            else
                rho_air_a=in{2};
                mu_air_a=in{3};
                cp_air_a=in{4};
                k_air_a=in{5};
            end
            
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
        
        function h=htube2htf_simple(HTF,T,MassFR,r_u,l)
            mu=HTF.mu(T);
            k=HTF.k(T);
            cp=HTF.cp(T);
            Re=2.*MassFR./(pi*r_u.*mu);
            Pr=cp.*mu./k;
            Gz=2*r_u./l.*Re.*Pr;
            Nu=zeros(size(Re));
            ilam=Re<3000;
            Nu(ilam)=4.36+0.0668*Gz(ilam)./(1+0.04*Gz(ilam).^(2/3));
            f=(0.79*log(Re(~ilam))-1.64).^(-2);
            Nu(~ilam)=(f/8.*(Re(~ilam)-1000).*Pr(~ilam))/(1+12.7*(f/8).^0.5.*(Pr(~ilam).^(2/3)-1));
            Nu(~ilam)=Nu(~ilam)*(1+1/(l/2/r_u)^(2/3));
            h=Nu.*k./r_u/2;  %W/m2-K
        end
        
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
            a=air();
            cp=a.cp(T_avg);
            mu=a.mu(T_avg);
            k=a.k(T_avg);
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
        
        function h=htube2htf(HTF,T,MassFR,r_u,l)
            mu=HTF.mu(T);
            k=HTF.k(T);
            cp=HTF.cp(T);
            Re=2.*MassFR./(pi*r_u.*mu);
            Pr=cp.*mu./k;
            Gz=2*r_u./l.*Re.*Pr;
            Nu=zeros(size(Re));
            ilam=Re<3000;
            Nu(ilam)=4.346*(1+(Gz(ilam)/29.6).^2).^(1/6).*(1+(Gz(ilam)/19.04./(1+(Pr(ilam)/0.0207).^(2/3)).^(1/2)./(1+(Gz(ilam)/29.6).^2).^(1/3)).^(3/2)).^(1/3);
            f=(0.79*log(Re(end))-1.64).^(-2);
            Nu_fd=(f/8.*(Re(end)-1000).*Pr(end))./(1+12.7*(f/8).^0.5.*(Pr(end).^(2/3)-1));
            % Nu_fd=mean(Nu(~ilam));
            dCDdx=(0.68+3000*Re(~ilam).^(-0.81)).*Pr(~ilam).^(-1/6)*0.1.*(l(~ilam)/2/r_u).^(-0.9);
            Nu(~ilam)=Nu_fd*(1+dCDdx);
            h=Nu.*k/r_u/2;  %W/m2-K
        end

    end
end