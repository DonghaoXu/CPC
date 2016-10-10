function [ref,inters,ia,oi]=refl3D(inc,stp,yl2,yl1,yr2,yr1,r_g,t_g,r_a,gap,tol,theta)
% ref is the vector showing reflected ray's direction
% inters is the intersection of reflected rays with the surface
% inc is the vector showing the incident ray's direction
% sp is the start point
% global yl2 yl1 yr2 yr1 yc ya s R ina r
%rotation matrix R=[cosa, -sina; sina cosa] clockwisely rotate by a (rad)
%y(1x2)=x(1x2)R
global ina s
ina=inc;
s=stp;
inc2D = inc(1: 2);
inc2D = normv(inc2D);
R=[inc2D(1),-inc2D(2);inc2D(2),inc2D(1)];
iR=[inc2D(1),inc2D(2);-inc2D(2),inc2D(1)];
ryl2=(yl2-ones(size(yl2,1),1)*stp(1: 2))*R;
ryl1=(yl1-ones(size(yl1,1),1)*stp(1: 2))*R;
ryr2=(yr2-ones(size(yr2,1),1)*stp(1: 2))*R;
ryr1=(yr1-ones(size(yr1,1),1)*stp(1: 2))*R;
% ryc=(yc-ones(size(yc,1),1)*s)*R;
% rya=(ya-ones(size(ya,1),1)*s)*R;
rs=-stp(1: 2)*R;

%ryl2=ryl2(ryl2(:,1)>0,:);
%ryl1=ryl1(ryl1(:,1)>0,:);
%ryr2=ryr2(ryr2(:,1)>0,:);
%ryr1=ryr1(ryr1(:,1)>0,:);
% ryc=ryc(ryc(:,1)>=0,:);
% rya=rya(rya(:,1)>0,:);

%get the point in the profile vector that is the closest to the
%intersection point
[il2,il2s]=getInter(ryl2);
[il1,il1s]=getInter(ryl1);
[ir2,ir2s]=getInter(ryr2);
[ir1,ir1s]=getInter(ryr1);
% [ia,ias]=getInter(rya,5,s,R,r);
ig=[inf,0];
igi=[inf,0];
ic=[inf,0];
if t_g~=0
    if abs(rs(2))<=r_g
        g1=rs(1)-sqrt(r_g^2-rs(2)^2);
        g2=rs(1)+sqrt(r_g^2-rs(2)^2);
        if g1>10e-8
            ig=[g1,0];
        elseif g2>10e-8
            ig=[g2,0];
        end
        if abs(rs(2))<=r_g-t_g
            gi1=rs(1)-sqrt((r_g-t_g)^2-rs(2)^2);
            gi2=rs(1)+sqrt((r_g-t_g)^2-rs(2)^2);
            if gi1>10e-8
                igi=[gi1,0];
            elseif gi2>10e-8
                igi=[gi2,0];
            end
            if abs(rs(2))<=r_a && rs(1)>0
                ic=[rs(1)-sqrt(r_a^2-rs(2)^2),0];
            end
        end
    end
else
	if abs(rs(2))<=r_a && rs(1)>0
        ic=[rs(1)-sqrt(r_a^2-rs(2)^2),0];
	end
end
igs=ig;
igis=igi;
ics=ic;
% [ic,ics]=getInter(ryc,6,s,R,r);

%find the shortest path and rotate it back
allinter=[il2,il2s;il1,il1s;ir1,ir1s;ir2,ir2s;inf,inf,inf,inf;ic,ics;ig,igs;igi,igis];
%ind_tooclose1=(allinter(:,1)<0.0035*t_g);
%ind_tooclose2=(allinter(:,3)<0.0035*t_g);
%ind_tooclose=ind_tooclose1 | ind_tooclose2;
%allinter(ind_tooclose,:)=inf;
%ind_tooclose1=(allinter(:,1)<0.08*r_g);
%ind_tooclose2=(allinter(:,3)<0.08*r_g);
%ind_tooclose=ind_tooclose1 | ind_tooclose2;
%ind23=([1;1;1;1;0;0;0;0]>0);
%ind_tooclose=ind_tooclose & ind23;
%allinter(ind_tooclose,:)=inf;
if any(any(isfinite(allinter(:,[1,3]))))
    [~,ind]=min(allinter(:,1));
    inters1=allinter(ind,1:2)*iR+stp(1: 2);
    inters2=allinter(ind,3:4)*iR+stp(1: 2);
    if abs(inters1(1)-inters2(1))<1e-5
        inters=(inters1+inters2)/2;
    else
%refine the intersection by solving ODE
        inters=refint(inters1,ind,theta,r_a,r_g,gap);
    end
    zdist = norm(stp(1: 2) - inters) * inc(3) / norm(inc(1: 2));
    inters = [inters, stp(3) + zdist];

%get the reflected direction
    [ref,ia,oi]=refdir3D(inc,inters,ind,theta,r_a);
else
    ref=[NaN,NaN,NaN];
    inters=[NaN,NaN,NaN];
    ia=NaN;
    oi=0;
end

end