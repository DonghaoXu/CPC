function [s,spr,spl]=sprefl(r_a,r_g,gap,theta)
angle_i=asin(r_a/(r_g+gap));
l_involute=sqrt((r_g+gap)^2-r_a^2)+(angle_i+theta)*r_a; %length of involute
spr=[r_a*cos(theta)+l_involute*sin(theta), r_a*sin(theta)-l_involute*cos(theta)];
spl=[-spr(1),spr(2)];
s=spl(1);
end