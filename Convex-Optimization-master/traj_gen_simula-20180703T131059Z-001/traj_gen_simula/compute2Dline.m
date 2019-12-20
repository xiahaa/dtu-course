function [a,b,c] = compute2Dline(x1, x2)
%% from inhomogeneous to homogeneous
x1h = [x1';1];
x2h = [x2';1];
l = cross(x1h,x2h);
a = l(1)*10;
b = l(2)*10;
c = l(3)*10;
end