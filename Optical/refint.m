function y=refint(inters1,flag,theta,r_a,r_g,gap)
switch(flag)
    case 1
%         opts=odeset('Events',@eaten_event,'AbsTol',tol);
%         [~,yl2]=ode45(@(t,y)profilel2(t,y,theta,r_a),[inters1(1),inters2(1)],inters1,opts);
%         y=yl2(end,:);
        y=Newton_DX(@(p) y2xy(p,theta,r_a,r_g,gap),inters1);
    case 2
%         opts=odeset('Events',@eaten_event,'AbsTol',tol);
%         [~,yl1]=ode45(@(t,y)profilel1(t,y,r_a),[inters1(1),inters2(1)],inters1,opts);
%         y=yl1(end,:);
%         options = optimoptions('fsolve','Jacobian','on','Display','off','DerivativeCheck','off');
%         y=fsolve(@(p) y1xy(p,r_a,r_g,gap),inters1,options);
        y=Newton_DX(@(p) y1xy(p,r_a,r_g,gap),inters1);
    case 3
%         opts=odeset('Events',@eaten_event,'AbsTol',tol);
%         [~,yr1]=ode45(@(t,y)profiler1(t,y,r_a),[inters1(1),inters2(1)],inters1,opts);
%         y=yr1(end,:);
%         options = optimoptions('fsolve','Jacobian','on','Display','off','DerivativeCheck','off');
%         y=fsolve(@(p) y1xy(p,r_a,r_g,gap),inters1,options);
        y=Newton_DX(@(p) y1xy(p,r_a,r_g,gap),inters1);
    case 4
%         opts=odeset('Events',@eaten_event,'AbsTol',tol);
%         [~,yr2]=ode45(@(t,y)profiler2(t,y,theta,r_a),[inters1(1),inters2(1)],inters1,opts);
%         y=yr2(end,:);
        y=Newton_DX(@(p) y2xy(p,theta,r_a,r_g,gap),inters1);
    otherwise
        y=inters1;
end
end

% function [val,isterm,dir]=eaten_event(t,y)
% global s ina
% dir=0;
% isterm=1;
% val=(y(2)-s(2))*ina(1)-(y(1)-s(1))*ina(2);
% end

function [f,J]=y1xy(p,r_a,r_g,gap)
global s ina
x=p(1); y=p(2);
d2=x^2+y^2;
f=zeros(2,1);
f(1)=(y-s(2))*ina(1)-(x-s(1))*ina(2);
beta=atan(-y/abs(x));
alpha=asin(r_a/sqrt(d2));
theta=pi/2-alpha-beta;
theta0=-asin(r_a/(r_g+gap));
l=sqrt((gap+r_g)^2-r_a^2)+(theta-theta0)*r_a;
f(2)=d2-r_a^2-l^2;
plpx=r_a*(r_a/cos(alpha)*x/d2^(3/2)-cos(beta)^2*sign(x)*y/x^2);
plpy=r_a*(r_a/cos(alpha)*y/d2^(3/2)+cos(beta)^2/abs(x));
J=[-ina(2),ina(1);2*x-2*l*plpx,2*y-2*l*plpy];
end

function [f,J]=y2xy(p,thetac,r_a,r_g,gap)
global s ina
x=p(1); y=p(2);
f=zeros(2,1);
f(1)=(y-s(2))*ina(1)-(x-s(1))*ina(2);
l2=sqrt((gap+r_g)^2-r_a^2)+(thetac+asin(r_a/(r_g+gap)))*r_a;
d2=x^2+y^2;
beta=atan(-y/abs(x));
alpha=asin(r_a/sqrt(d2));
theta=pi/2-alpha-beta;
OM=r_a/sin(thetac);
NS=(y-OM)*tan(thetac);
MS=NS/sin(thetac);
TM=r_a*cot(thetac);
AS=l2+TM+MS;
BS=(NS+abs(x))*sin(thetac);
AB=AS-BS;
lp=(theta-thetac)*r_a+l2+AB;
f(2)=d2-r_a^2-lp^2;
plpx=r_a*(r_a/cos(alpha)*x/d2^(3/2)-cos(beta)^2*sign(x)*y/x^2)-sign(x)*sin(thetac);
plpy=r_a*(r_a/cos(alpha)*y/d2^(3/2)+cos(beta)^2/abs(x))+cos(thetac);
J=[-ina(2),ina(1);2*x-2*lp*plpx,2*y-2*lp*plpy];
end

function x=Newton_DX(fh,x0)
tol=1e-8;
x=x0';
xold=x;
while 1
    [f,J]=feval(fh,x);
    if norm(f)<tol, break; end
    x=xold-J\f;
    rx=(x-xold)./(xold+eps);
    if norm(rx)<tol, break; end
    xold=x;
end
x=x';
end