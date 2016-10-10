function [ref,ia,oi]=refdir(inc,y,flag,theta,r)
%by giving incident angle and intersection and flag showing which curve it
%is, it returns the reflected direction which is a unit vector
% global theta r
ri=1.5167;
% 
% function y=fr(x)
% if x==1
%     y=1;
% else
%     y=ri;
% end
% end

switch(flag)
    case 1
        %slope at the intersection
        t=sqrt(y(1)^2+y(2)^2-r^2);
        sina=(-t*y(2)-r*y(1))/(y(1)^2+y(2)^2);
        cosa=(r*y(2)-t*y(1))/(y(1)^2+y(2)^2);
        k=[(cos(theta)+sina),-(sin(theta)+cosa)];
        k=k/norm(k);
        %reflected direction
        ref=2*(k*inc')*k-inc;
        ia=asin(abs(k*inc'))/pi*180;
        oi=1;
    case 2
        %slope at the intersection
        sina=(r*y(2)-sqrt(y(1)^2+y(2)^2-r^2)*y(1))/(y(1)^2+y(2)^2);
        cosa=(r*y(1)+sqrt(y(1)^2+y(2)^2-r^2)*y(2))/(y(1)^2+y(2)^2);
        k=[cosa,sina];
        k=k/norm(k);
        %reflected direction
        ref=2*(k*inc')*k-inc;
        ia=asin(abs(k*inc'))/pi*180;
        oi=2;
    case 3
        %slope at the intersection
        sina=(sqrt(y(1)^2+y(2)^2-r^2)*y(1)+r*y(2))/(y(1)^2+y(2)^2);
        cosa=(r*y(1)-sqrt(y(1)^2+y(2)^2-r^2)*y(2))/(y(1)^2+y(2)^2);
        k=[cosa,sina];
        k=k/norm(k);
        %reflected direction
        ref=2*(k*inc')*k-inc;
        ia=asin(abs(k*inc'))/pi*180;
        oi=3;
    case 4
        %slope at the intersection
        t=sqrt(y(1)^2+y(2)^2-r^2);
        sina=(t*y(2)-r*y(1))/(y(1)^2+y(2)^2);
        cosa=(r*y(2)+t*y(1))/(y(1)^2+y(2)^2);
        k=[cos(theta)-sina,sin(theta)+cosa];
        k=k/norm(k);
        %reflected direction
        ref=2*(k*inc')*k-inc;
        ia=asin(abs(k*inc'))/pi*180;
        oi=4;
    case 6
        %slope at the intersection
        n=[y(1),y(2)];
        n=n/norm(n);
        %reflected direction
        ref=inc-2*(n*inc')*n;
        ia=acos(abs(n*inc'))*180/pi;
        oi=-1;
    case 5
        ref=[NaN,NaN];
        oi=0;
    case 7  %refraction at outer surface of glass
        n=[y(1),y(2)];
        n=n/norm(n);
        cosia=max(-1,min(1,inc*n'));
        sinia=sqrt(1-cosia^2);
        if cosia<0 %into glass
            sinra=sinia/ri;
            cosra=sqrt(1-sinra^2);
%             ref=inc-(cosra*fr(cosia)+cosia)*n;
            ref=inc-(cosra*ri+cosia)*n;
            oi=7;
        else %out of glass
            sinra=sinia*ri;
            if sinra<1
                cosra=sqrt(1-sinra^2);
%                 ref=inc-(cosia-cosra/fr(cosia))*n;
                ref=inc-(cosia-cosra/ri)*n;
            else
                ref=inc-2*(n*inc')*n;
            end
            oi=-7;
        end
        ref=ref/norm(ref);
        ia=acos(abs(cosia))*180/pi;
    case 8 %refraction at inner surface of glass
        n=[y(1),y(2)];
        n=n/norm(n);
        cosia=max(-1,min(1,inc*n'));
        sinia=sqrt(1-cosia^2);
        if cosia>0 %into glass
            sinra=sinia/ri;
            cosra=sqrt(1-sinra^2);
%             ref=inc+(cosra*fr(cosia)-cosia)*n;
            ref=inc+(cosra*ri-cosia)*n;
            oi=8;
        else %out of glass
            sinra=sinia*ri;
            if sinra<1
                cosra=sqrt(1-sinra^2);
%                 ref=inc-(cosia+cosra/fr(cosia))*n;
                ref=inc-(cosia+cosra/ri)*n;
                oi=-8;
            else     %total internal reflection
                ref=inc-2*(n*inc')*n;
                oi=8*sqrt(-1);
            end
        end
        ref=ref/norm(ref);
        ia=acos(abs(cosia))*180/pi;
end
        
end

