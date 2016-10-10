function [y1,y2,alpha,phi,delta,h,AST]=IncAngle(n,Time,LAT,LONG,LSM,beta,gama,varargin)
%n is the nth day of the year
%Time is local standard time during DST
%LAT is local latitude
%LONG is local longitude
%LSM is local standard meridians
%beta is tilted angle of the surface
%gama is the surface azimuth
%varargin: empty is no DST correction; any is yes
n=n(:); Time=Time(:);
%Equation of time
ET=(9.87*sin(4*pi*(n-81)/364)-7.53*cos(2*pi*(n-81)/364)-1.5*sin(2*pi*(n-1)/364))/60;

AST=Time-(~isempty(varargin))+ET+(LONG-LSM)/15;    %correction for DST

if isempty(beta) || isempty(gama)
    y1=[];y2=[];alpha=[];phi=[];delta=[];h=[];
    return
end

%solar declination
delta=23.45*sin(2*pi*(284+n)/360)*pi/180;

%hour angle
h=15*(AST-12)*pi/180;

%solar altitude
alpha=asin(cos(LAT*pi/180).*cos(delta).*cos(h)+sin(LAT*pi/180).*sin(delta));

%solar azimuth
phi=acos((sin(alpha).*sin(LAT*pi/180)-sin(delta))./cos(alpha)./cos(LAT*pi/180)).*sign(h);
% In the morning, phi<0; in the afternoon, phi>0

%incident angle
y2=acos(cos(alpha).*cos(phi-gama*pi/180).*sin(beta*pi/180)+sin(alpha).*cos(beta*pi/180));

y1=y2*180/pi;
alpha=alpha*180/pi;
phi=phi*180/pi;
delta=delta*180/pi;
h=h*180/pi;
end