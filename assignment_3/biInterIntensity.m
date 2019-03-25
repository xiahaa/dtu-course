
function Icurve = biInterIntensity(im, curve)
    point1 = floor(curve);point2 = point1;point3 = point1;
    point2(2,:) = point2(2,:) + 1;
    point3(1,:) = point3(1,:) + 1;
    point4 = point2;
    point4(1,:) = point4(1,:) + 1;
    
    d1 = curve(1,:) - point1(1,:);
    d2 = curve(2,:) - point1(2,:);
    
    s1 = (1-d1).*(1-d2);
    s2 = (1-d1).*d2;
    s3 = d1.*(1-d2);
    s4 = (d1).*(d2);
    
    index1 = sub2ind(size(im), point1(1,:),point1(2,:));
    index2 = sub2ind(size(im), point2(1,:),point2(2,:));
    index3 = sub2ind(size(im), point3(1,:),point3(2,:));
    index4 = sub2ind(size(im), point4(1,:),point4(2,:));
    
    Icurve = im(index1).*s1 + ...
             im(index2).*s2 + ...
             im(index3).*s3 + ...
             im(index4).*s4;
end