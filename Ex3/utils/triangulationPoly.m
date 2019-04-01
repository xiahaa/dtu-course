function P = triangulationPoly(x1,P1,x2,P2)
% TODO check this function, may be not correct
    % F
    M1 = P1(1:3,1:3);
    c1 = -M1\P1(1:3,4);
    
    M2 = P2(1:3,1:3);
    c2 = -M2\P2(1:3,4);
    
    e1 = P2 * [c1;1];
    e2 = P1 * [c2;1];
    
    P1inv = pinv(P1);
    
    se = [0 -e1(3) e1(2);e1(3) 0 -e1(1);-e1(2) e1(1) 0];
    F = se*P2*P1inv;

    % T
    theta1 = atan2(-(e1(2)-e1(3)*x1(2)),e1(1)-e1(3)*x1(1));
    theta2 = atan2(-(e2(2)-e2(3)*x2(2)),e2(1)-e2(3)*x2(1));
    
    T1 = [cos(theta1),-sin(theta1),0;sin(theta1),cos(theta1),0;0,0,1] * [1 0 -x1(1);0 1 -x1(2);0 0 1];
    T2 = [cos(theta2),-sin(theta2),0;sin(theta2),cos(theta2),0;0,0,1] * [1 0 -x2(1);0 1 -x2(2);0 0 1];
    
    Fn = inv(T2)'*F*inv(T1);
    
    a = Fn(2,2); b = Fn(2,3); c = Fn(3,2); d = Fn(3,3);
    f1 = -Fn(2,1)/b; f2 = -Fn(1,2)/c;
    
    %cs = [ t^6, t^5, t^4, t^3, t^2, t, 1]
    poly = [-a*c*f1^4*(a*d - b*c), (a^2 + c^2*f2^2)^2 - a*d*f1^4*(a*d - b*c) - b*c*f1^4*(a*d - b*c), ...
            2*(a^2 + c^2*f2^2)*(2*c*d*f2^2 + 2*a*b) - 2*a*c*f1^2*(a*d - b*c) - b*d*f1^4*(a*d - b*c), ...
            2*(a^2 + c^2*f2^2)*(b^2 + d^2*f2^2) + (2*c*d*f2^2 + 2*a*b)^2 - 2*a*d*f1^2*(a*d - b*c) - 2*b*c*f1^2*(a*d - b*c), ...
            2*(b^2 + d^2*f2^2)*(2*c*d*f2^2 + 2*a*b) - a*c*(a*d - b*c) - 2*b*d*f1^2*(a*d - b*c), ...
            (b^2 + d^2*f2^2)^2 - a*d*(a*d - b*c) - b*c*(a*d - b*c), -b*d*(a*d - b*c)];
    
    rs = roots(poly);
    
    mins = 1e6;
    t = 0;
    for i = 1:size(rs,1)
        if abs(imag(rs(i)))<1e-6
            if evals(real(rs(i)),a,b,c,d,f1,f2) < mins
                t = real(rs(i));
            end
        end
    end
 
    % 
    l1 = [t*f1;1;-t];
    l2 = [-f2*(c*t+d);a*t+b;c*t+d];
    
    l3 = [-l1(2);l1(1);0];
    l4 = [-l2(2);l2(1);0];
    
    u1 = cross(l1,l3);u1 = u1./u1(3);
    u2 = cross(l2,l4);u2 = u2./u2(3);
    
    xu1 = inv(T1)*u1;
    xu2 = inv(T2)*u2;
    
    P = triangulationEigen(xu1,P1,xu2,P2);
end

function s = evals(t,a,b,c,d,f1,f2)
    s = t^2/(1+(t*f1)^2) + (c*t+d)^2/((a*t+b)^2+f2^2*(c*t+d)^2);
end