
function P = triangulationMidpoint(x1,P1,x2,P2)
%Implementation of mid-point triangulation method.
    M1 = P1(1:3,1:3);
    c1 = -M1\P1(1:3,4);
    
    M2 = P2(1:3,1:3);
    c2 = -M2\P2(1:3,4);
    
    % ray
    r1 = M1\x1;
    r2 = M2\x2;
    
    A = [r1 -r2];
    b = c2 - c1;
    
    x = A\b;
    
    P = (c1 + x(1).*r1 + c2 + x(2).*r2)*0.5;
end
