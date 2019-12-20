clc;close all;clear all;

[x,y,g] = gaussian(2, 5);
surf(x,y,g);
[x1,y1,g1] = dgaussian(2, 5, 1);
[x2,y2,g2] = dgaussian(2, 5, 2);
surf(x1,y1,g1);
surf(x2,y2,g2);


function [x,y,g] = gaussian(t, sigma)
    m = round(t*sigma);
    [X,Y] = meshgrid(-m:1:m,-m:1:m);
    f = @(x,y) (1/(2*pi*t*t).*exp(-0.5/(t*t).*(x.*x+y.*y)));
    g = f(X,Y);
    x = -m:1:m;
    y = -m:1:m;
end

function [x,y,g] = dgaussian(t, sigma, dir)
    m = round(t*sigma);
    [X,Y] = meshgrid(-m:1:m,-m:1:m);
    dfx = @(x,y) (1/(2*pi*t*t).*exp(-0.5/(t*t).*(x.*x+y.*y)).*(-1/(t*t).*x));
    dfy = @(x,y) (1/(2*pi*t*t).*exp(-0.5/(t*t).*(x.*x+y.*y)).*(-1/(t*t).*y));
    if dir == 1
       g = dfx(X,Y);
    else
       g = dfy(X,Y);
    end
    x = -m:1:m;
    y = -m:1:m;
end