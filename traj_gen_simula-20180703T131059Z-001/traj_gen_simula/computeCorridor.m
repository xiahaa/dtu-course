function [xl,yl,xr,yr] = computeCorridor(a,b,c,x,y,l)

% if abs(b) < 1e-6
%     xl = x - l;
%     xr = x + l;
%     yl = y;
%     yr = y;
% else
%     cl = c - l*sqrt(a*a+b*b);
%     cr = c + l*sqrt(a*a+b*b);
%     xl = x;
%     xr = x;
%     yl = -1/b*(a*x+cl);
%     yr = -1/b*(a*x+cr);
% end

%% orthogonal direction
v1 = [a;b];
v1 = v1./norm(v1);
p1 = [x;y]+l.*v1;
p2 = [x;y]-l.*v1;
xl = p1(1);
yl = p1(2);
xr = p2(1);
yr = p2(2);

end