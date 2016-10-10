function y=RayDisp(filename,t,hd)
%t is the angle vector user wants to display
%hd=1 wait for keyboard to continue
%hd=0 display all
%hd=other dont delete
if ischar(filename)
    load(filename,'history')
else
    history=filename;
end
% tol=1e-5;
r_a=history.RadReceiver;
r_g=history.RadGlass;
t_g=history.ThickGlass;
w=history.Width;
gap=history.Gap;
theta_deg=history.AccHalfAng;
theta=theta_deg*pi/180;
[~,~,spl]=sprefl(r_a,r_g,gap,theta);
% [yl1,yl2,yr1,yr2,yc,yg_o,yg_i,~]=profileXCPC(theta,r_a,w,r_g,t_g,gap,spl,tol);
[yl1,yl2,yr1,yr2,yc,yg_o,yg_i,~]=profileAnaXCPC(theta,r_a,w,r_g,t_g,gap,spl);
plot(yl1(:,1),yl1(:,2),'color',[38 188 213]/255)
hold on
axis equal
plot(yl2(:,1),yl2(:,2),'color',[38 188 213]/255)
plot(yr1(:,1),yr1(:,2),'color',[38 188 213]/255)
plot(yr2(:,1),yr2(:,2),'color',[38 188 213]/255)
plot(yc(:,1),yc(:,2),'color',[254 67 101]/255)
plot(yg_o(:,1),yg_o(:,2),'color',[29 191 151]/255)
plot(yg_i(:,1),yg_i(:,2),'color',[29 191 151]/255)
xl=xlim;
yl=ylim;
xlim(xl*1.1);
ylim(yl*1.1);
set(gcf,'color',[1 1 1])
if nargin==2 && any(isnan(t))
    return
end
[l,~,na]=size(history.x);
h=zeros(size(1:50:l));
if nargin<2
    t=1:na;
end
step=30;
for j=1:numel(t)
    s=((t(j)>0)-0.5)*2;
    for i=1:step:l
        sd=s*history.x(i,~isnan(history.x(i,:,abs(t(j))+1)),abs(t(j))+1);
        sy=history.y(i,~isnan(history.y(i,:,abs(t(j))+1)),abs(t(j))+1);
        h((i-1)/step+1)=plot(sd,sy,'k');
        hold on
    end
    title(sprintf('Incidence angle = %.1f Deg.', t(j)))
    pause(0.01)
    if hd==1
        keyboard
        delete(h)
    elseif hd==0
        delete(h)
    end  
end
y=gcf;
end
