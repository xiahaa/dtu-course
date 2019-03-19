function I = drawline(I,x,y,color)
    M = size(I,1);N = size(I,2);
    nPoints = max(abs(diff(x)), abs(diff(y)))+1;  % Number of points in line
    rIndex = round(linspace(y(1), y(2), nPoints));  % Row indices
    cIndex = round(linspace(x(1), x(2), nPoints));  % Column indices
    index = sub2ind([M,N], rIndex, cIndex);     % Linear indices
    nn = [-M +M -1 +1];
    index = [index index+nn(1) index+nn(2) index+nn(3) index+nn(4)];
    index = index(index > 0 & index < M*N);
    I(index) = color(1);  % Set the line pixels to the max value of 255 for uint8 types
    I(index+M*N) = color(2);
    I(index+M*N*2) = color(3);
end