function t = r2t(r, n)
if n < 3
    t = nan;
    return
end
t = r.*sqrt((n - 2)./(1 - r.^2));
end