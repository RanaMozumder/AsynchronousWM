function [T_vec, T, theta] = CalculateT(SpkCount,angles)
    num = 0;
    den = 0;
    SpkCount(find(isnan(SpkCount)))=0;
    if isempty(angles)
        angles = [0 45 90 135 180 225 270 315];
    end
    for class=1:length(SpkCount)
        num = num + SpkCount(class)*exp(1i * deg2rad(angles(class)));
        den = den + SpkCount(class);
    end
    T_vec = num/den;
    T = abs(T_vec);
    theta = rad2deg(angle(T_vec));
    if theta<0
        theta = theta + 360;
    end
%     if isnan(theta)
%         pause
%     end
end