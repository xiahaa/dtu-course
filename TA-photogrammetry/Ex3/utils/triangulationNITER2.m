
function P = triangulationNITER2(x1,P1,x2,P2)
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
    
    S = [1 0 0;0 1 0];
    Ftilde = S*F*S';
    
    n2 = S*F*x1;
    n1 = S*F'*x2;
    
    a = n2'*Ftilde*n1;
    b = 0.5*(n2'*n2+n1'*n1);
    c = x2'*F*x1;
    d = sqrt(b*b-a*c);
    lamda = c/(b+d);
    
    dx2 = lamda * n2;
    dx1 = lamda * n1;
    
    n2 = n2 - Ftilde * dx1;
    n1 = n1 - Ftilde' * dx2;
    
    lamda = lamda * 2*d / (n2'*n2+n1'*n1);
    dx2 = lamda * n2;
    dx1 = lamda * n1;
    
    xu2 = x2 - S'*dx2;
    xu1 = x1 - S'*dx1;
    
    P = triangulationEigen(xu1,P1,xu2,P2);
end