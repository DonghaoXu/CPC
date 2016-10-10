function y=IncAng2D(alpha,phi,beta)
% alpha: Solar altitude
% phi: Solar azimuth
alpha=alpha(:)/180*pi;
phi=phi(:)/180*pi;
beta=beta/180*pi;
% gama=gama/180*pi;
p=[cos(-phi)';sin(-phi)';tan(alpha)'];
R=[cos(beta),0,-sin(beta);0,1,0;sin(beta),0,cos(beta)];
p=R*p;
cosia2=[0,1]*p(2:3,:)./sqrt(p(2,:).^2+p(3,:).^2);
y=acos(cosia2')*180/pi.*sign(phi);
end