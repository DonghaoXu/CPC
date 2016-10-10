function h=htube2htf(T,MassFR,r_u,l)
mu=mu_HTF(T);
k=k_HTF(T);
cp=cp_HTF(T);
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

% nu=4.346*(1+(gz/29.6).^2).^(1/6).*(1+(gz/19.04./(1+(Pr/0.0207).^(2/3)).^(1/2)./(1+(gz/29.6).^2).^(1/3)).^(3/2)).^(1/3)