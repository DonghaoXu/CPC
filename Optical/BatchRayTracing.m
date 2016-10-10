function [history,time]=BatchRayTracing(theta_deg,r_a,r_g,t_g,w,gap,iai,l)
%global theta r w yl2 yl1 yr2 yr1 yc ya spr spl ipl ipr tol
tol=1e-5;
%r_a=2.9;  %radius of absorber: cm
%r_g=3.25; %radius of glass envelope: cm
%t_g=0.1;  %thickness of glass envelope: cm
%w=21.4; %width of aperture:cm
%h=11 height of XCPC: cm for validation
%gap=0.4;  %gap between glass envelope and XCPC: cm
%r_ut=0.6; %radius of u-tube
% r_g=r_g1/r_a1;
% t_g=t_g1/r_a1;
% w=w1/r_a1;
% gap=gap/r_a1;
% r_a=1;
theta=pi/180*theta_deg;

%%
%starting points
[~,~,spl]=sprefl(r_a,r_g,gap,theta);

%%
%profiles
[yl1,yl2,yr1,yr2,~,~,~,~]=profileAnaXCPC(theta,r_a,w,r_g,t_g,gap,spl);
%[yl1,yl2,yr1,yr2,~,~,~,~]=profileXCPC(theta,r_a,w,r_g,t_g,gap,spl,tol);

%%
% GPUarray
% yl2=gpuArray(yl2); yl1=gpuArray(yl1); yr2=gpuArray(yr2); yr1=gpuArray(yr1);
% r_g=gpuArray(r_g); t_g=gpuArray(t_g); r_a=gpuArray(r_a); gap=gpuArray(gap);
% l=gpuArray(l); iai=gpuArray(iai); tol=gpuArray(tol); theta=gpuArray(theta);

%%
[history,time]=RayTracing(yl2,yl1,yr2,yr1,r_g,t_g,r_a,gap,l,iai,tol,theta);

%%
% history.x=history.x*r_a1;
% history.y=history.y*r_a1;
history.Width=w;
history.RadReceiver=r_a;
history.RadGlass=r_g;
history.ThickGlass=t_g;
history.Gap=gap;
history.CR=w/2/pi/r_a;
history.AccHalfAng=theta_deg;

end