function pos = bancroft_fast(pr, sat_pos, xinit)
    x0 = xinit(1);
    y0 = xinit(2);
    z0 = xinit(3);
    cdt = 299792458 * xinit(4);

    x1 = sat_pos(:,1);
    y1 = sat_pos(:,2);
    z1 = sat_pos(:,3);
        
    B = [x1 y1 z1 pr];
    
    a = 0.5.*(x1.^2+y1.^2+z1.^2-pr.^2);
    e = ones(size(pr,1),1);
    
    if size(pr,1) == 4
        BBB = inv(B);
    else
        BBB = inv(B'*B)*B';
    end
    
    v1 = BBB * e;
    v2 = BBB * a;
        
    a1 = v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3)-v1(4)*v1(4);
    a2 = v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)-v1(4)*v2(4);
    a3 = v2(1)*v2(1)+v2(2)*v2(2)+v2(3)*v2(3)-v2(4)*v2(4);
    
    sqrtm = sqrt((a2-1)^2-a1*a3);
    
    z1 = (-(a2-1)+sqrtm)/a1;
    z2 = (-(a2-1)-sqrtm)/a1;
    
    M = eye(4);
    M(4,4) = -1;
    x01 = M*BBB*(z1.*e+a);
    x02 = M*BBB*(z2.*e+a);
    
    dist1 = norm(x01(1:3),2);
    dist2 = norm(x02(1:3),2);
    if dist1 > dist2
        pos = x02;
    else
        pos = x01;
    end
end