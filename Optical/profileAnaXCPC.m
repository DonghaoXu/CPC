function [yl1,yl2,yr1,yr2,yc,yg_o,yg_i,ya]=profileAnaXCPC(theta,r_a,w,r_g,t_g,gap,spl,varargin)
theta0=-asin(r_a/(r_g+gap));
if isempty(varargin)
    options=optimoptions('fsolve','Display','off');
else
    options=varargin{1};
end
if -w/2<spl(1)
    ta=[theta,theta-1e-5:-pi/180:theta0+1e-4]';
    yl1=zeros(size(ta,1),2);
    l=(ta-theta0)*r_a+sqrt((r_g+gap)^2-r_a^2);
    yl1(:,1)=-r_a*cos(ta)-l.*sin(ta);
    yl1(:,2)=r_a*sin(ta)-l.*cos(ta);
    sp = theta;
    while 1
        thetam=fsolve(@(ta) yl2xy(ta,r_a,r_g,gap,theta)+w/2,sp,options);
        ta=[theta:pi/180:thetam-1e-5,thetam]';
        yl2=zeros(size(ta,1),2);
        [yl2(:,1),yl2(:,2)]=yl2xy(ta,r_a,r_g,gap,theta);
        tyl2x = yl2(:, 1) + w/2 + 1e-5;
        ind = find(tyl2x < 0, 1);
        if isempty(ind)
            break
        else
            p = yl2(ind, :);
            sp = asin(p(2) / norm(p)) + acos(r_a / norm(p));
        end
    end
else
    theta_n=fsolve(@(theta) sprefl(r_a,r_g,gap,theta)+w/2,theta,options);
    yl2=[NaN,NaN];
    ta=[theta_n,theta_n-1e-5:-pi/180:theta0+1e-4]';
    yl1=zeros(size(ta,1),2);
    l=(ta-theta0)*r_a+sqrt((r_g+gap)^2-r_a^2);
    yl1(:,1)=-r_a*cos(ta)-l.*sin(ta);
    yl1(:,2)=r_a*sin(ta)-l.*cos(ta);    
end
yr2=[-yl2(:,1),yl2(:,2)];   %value=4
yr1=[-yl1(:,1),yl1(:,2)];   %value=3
agridc=(-pi:pi/50:pi)';    %angle grid for the tubular absorber
yc=[r_a*cos(agridc),r_a*sin(agridc)];   %value=6 for tubular absorber
yg_o=[r_g*cos(agridc),r_g*sin(agridc)];   %value=7 for glass envelope outside
yg_i=[(r_g-t_g)*cos(agridc),(r_g-t_g)*sin(agridc)];   %value=8 for glass envelope inside
ya=[yl2(end,:);(yl2(end,1)-r_g)/2,yl2(end,2);0,yl2(end,2);(yr2(end,1)+r_g)/2,yl2(end,2);yr2(end,:)]; %value=5 for aperture    
end

function [x,y]=yl2xy(ta,r_a,r_g,gap,theta)
cr1=(sin(ta)-sin(theta))/cos(theta);
cr21=cr1*sin(theta)-cos(theta);
cr22=cos(ta);
cr31=(cr21+cr22)*sin(theta);
cr32=-(pi/2+ta);
cr33=-cr1;
cr=cr31+cr32+cr33;
cl1=cos(ta)/cos(theta);
cl21=-cl1*sin(theta);
cl22=sin(ta);
cl31=(cl21+cl22)*sin(theta);
cl32=cl1;
cl=cl31+cl32+1;
g=r_g+gap-r_a;
l2=sqrt((r_a+g)^2-r_a^2)+r_a*(theta+asin(r_a/(r_a+g)));
S2=(pi/2+theta)*r_a;
l=(2*l2-S2-cr*r_a)./cl;
x=-r_a*cos(ta)-l.*sin(ta);
y=r_a*sin(ta)-l.*cos(ta);
%y=(Sx-lx)-(S2-l2)-RHS;
end