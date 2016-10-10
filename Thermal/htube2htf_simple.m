function h=htube2htf_simple(T,MassFR,r_u,l)
mu=mu_HTF(T);
k=k_HTF(T);
cp=cp_HTF(T);
Re=2.*MassFR./(pi*r_u.*mu);
Pr=cp.*mu./k;
Gz=2*r_u./l.*Re.*Pr;
Nu=zeros(size(Re));
ilam=Re<3000;
Nu(ilam)=4.36+0.0668*Gz(ilam)./(1+0.04*Gz(ilam).^(2/3));
f=(0.79*log(Re(~ilam))-1.64).^(-2);
Nu(~ilam)=(f/8.*(Re(~ilam)-1000).*Pr(~ilam))/(1+12.7*(f/8).^0.5.*(Pr(~ilam).^(2/3)-1));
Nu(~ilam)=Nu(~ilam)*(1+1/(l/2/r_u)^(2/3));
h=Nu.*k/r_u/2;  %W/m2-K
end

% nu=4.346*(1+(gz/29.6).^2).^(1/6).*(1+(gz/19.04./(1+(Pr/0.0207).^(2/3)).^(1/2)./(1+(gz/29.6).^2).^(1/3)).^(3/2)).^(1/3)