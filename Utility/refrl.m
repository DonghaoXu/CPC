function [rr, rl, iai, iarr, iarl] = refrl(i, n, nr)
i = normv(i(:)');
n = normv(n(:)');
cosia = i * n';
if cosia < 0
    n = -n;
    cosia = -cosia;
end
sinia = sqrt(1 - cosia^2);
iai = rad2deg(asin(sinia));

sinra = sinia / nr;
if sinra < 1
    cosra = sqrt(1 - sinra^2);
    rr = i + (nr * cosra - cosia) * n;
    iarr = rad2deg(asin(sinra));
    rr = normv(rr);
else
    rr = nan(size(i));
    iarr = nan;
end

rl = i - 2 * cosia * n;
iarl = iai;
rl = normv(rl);
end