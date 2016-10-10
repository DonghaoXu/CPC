function history=RayTracing_big(yl2,yl1,yr2,yr1,r_g,t_g,r_a,gap,l,iai,tol,theta)
assert(l == numel(iai))
iai = deg2rad(iai(:));
nray=10;
top = yl1(1, :);
if ~isnan(yl(1, 1))
    top = yl2(end, :);
end
pleft = top;
L = 100;
hit = zeros(l, 1);
sp = zeros(l, 2);

parfor i = 1: l
    lambda = rand();
    ia = iai(i);
    v1 = r_a * [cos(ia), sin(ia)];
    % v(2) = top(2). So ll should be as below
    ll = (top(2) - r_a * sin(ia)) / sin(-pi / 2 + ia);
    v2 = ll * [cos(-pi / 2 + ia); sin(-pi / 2 + ia)];
    v = v1 + v2;
    if v(1) > -top(1)
        pright = v;
    else
        pright = [-top(1), top(2)];
    end
    spi = pleft * lambda + pright * (1 - lambda) - v2 / ll * L;
    sp(i, :) = spi;
    fprintf('Incident angle: %.2f, Coordinates: [%.2f, %.2f]\n',...
        ia, spi(1), spi(2))
    % Start ray tracing
    k = 1;
    vra= v2 / ll;
    in = sp(i, :);
    while k<=nray
        [vra, in, ~, oi]=refl(vra,in,yl2,yl1,yr2,yr1,r_g,t_g,r_a,gap,tol,theta);
        if oi == -1
            hit(i) = 1;
            break
        elseif oi == 0
            break
        end
        k = k + 1;
    end
end

history.hit = hit;
history.sp = sp;
history.ia = iai;

end