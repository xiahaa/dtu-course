function Iw = polar_unwarping(I,cc)
    [M,N] = size(I);
    r1 = min(cc(1),M-cc(1));
    r2 = min(cc(2),N-cc(2));
    r = min(r1,r2);
    
    M1 = round(r);
    N1 = 360;
    
    rs = linspace(0,r,M1);
    as = linspace(0,2*pi,N1);
    
    Iw = zeros(M1,N1);
    
    for i = 1:numel(rs)
        rr = rs(i);
        for j = 1:numel(as)
            u = rr*cos(as(j)) + cc(2);
            v = rr*sin(as(j)) + cc(1);
            % interpolation
            u1 = (floor(u));
            v1 = (floor(v));
            % nearest neighbor interpolation
            if u1>0 && u1 < N && v1 >0 && v1 < M
                Iw(i,j) = I(v1,u1);
            end
        end
    end
    
    Iw = flipud(Iw);
end