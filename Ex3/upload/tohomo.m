function xh = tohomo(x)
%From inhomogeneous coordinate to homogeneous coordiante.
    xh = [x;ones(1,size(x,2))];
end