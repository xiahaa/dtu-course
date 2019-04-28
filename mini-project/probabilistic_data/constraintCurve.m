
function curve = constraintCurve(curve, m, n)
%Constrain curves being inside the image.
    id1 = curve(1,:) < 2;
    id2 = curve(1,:) > m-1;
    id3 = curve(2,:) < 2;
    id4 = curve(2,:) > n-1;
    curve(1,id1) = 2;
    curve(1,id2) = m - 1;
    curve(2,id3) = 2;
    curve(2,id4) = n - 1;
end