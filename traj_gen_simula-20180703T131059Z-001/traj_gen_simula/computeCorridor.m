function [xl,yl,xr,yr] = computeCorridor(a,b,c,x,y,l)

if abs(b) < 1e-6
    xl = x - l;
    xr = x + l;
    yl = y;
    yr = y;
else
    cl = c - l*sqrt(a*a+b*b);
    cr = c + l*sqrt(a*a+b*b);
    xl = x;
    xr = x;
    yl = -1/b*(a*x+cl);
    yr = -1/b*(a*x+cr);
end

% cl = c - l*sqrt(a*a+b*b);
% c2 = c + l*sqrt(a*a+b*b);
% %% l1,l2
% l1 = [a;b;c1];
% l2 = [a;b;c2];
%% norm

end