
function P = triangulationKNS(x1,P1,x2,P2)
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
    u = vec(F');

    E0 = 1e6;
    f0 = 1;% need this?
    
    xhat1 = x1(1); yhat1 = x1(2); 
    xhat2 = x2(1); yhat2 = x2(2);
    
    xtilde1 = 0; ytilde1 = 0;
    xtilde2 = 0; ytilde2 = 0;
    
    E = 0;
    iter = 1;
    maxIter = 50;
    while iter < maxIter
        
        xihat = [xhat1*xhat2 + xhat2*xtilde1 + xhat1*xtilde2; ...
                 xhat1*yhat2 + yhat2*xtilde1 + xhat1*ytilde2; ...
                 f0*(xhat1+xtilde1); ...
                 yhat1*xhat2 + xhat2*ytilde1 + yhat1*xtilde2; ...
                 yhat1*yhat2 + yhat2*ytilde1 + yhat1*ytilde2; ...
                 f0*(yhat1+ytilde1); ...
                 f0*(xhat2+xtilde2); ...
                 f0*(yhat2+ytilde2); ...
                 f0*f0];
        f02 = f0*f0;
        V0xi = [(xihat(3)^2+xihat(7)^2)/(f02) xihat(7)*xihat(8)/f02 xihat(7) xihat(3)*xihat(6)/f02 0 0 xihat(3) 0 0; ...
                xihat(7)*xihat(8)/f02 (xihat(3)^2+xihat(8)^2)/(f02) xihat(8) 0 xihat(3)*xihat(6)/f02 0 0 xihat(3) 0; ...
                xihat(7) xihat(8) xihat(9) 0 0 0 0 0 0; ...
                xihat(3)*xihat(6)/f02 0 0 (xihat(6)^2+xihat(7)^2)/(f02) xihat(7)*xihat(8)/f02 xihat(7) xihat(6) 0 0; ...
                0 xihat(3)*xihat(6)/f02 0 xihat(7)*xihat(8)/f02 (xihat(6)^2+xihat(8)^2)/(f02) xihat(8) 0 xihat(6) 0; ...
                0 0 0 xihat(7) xihat(8) xihat(9) 0 0 0; ...
                xihat(3) 0 0 xihat(6) 0 0 xihat(9) 0 0; ...
                0 xihat(3) 0 0 xihat(6) 0 0 xihat(9) 0; ...
                0 0 0 0 0 0 0 0 0];
    
        % update
        dot1 = dot(u,xihat);
        dot2 = dot(u,V0xi*u);
        
        xtilde1 = dot1 / dot2 * (u(1:3)'*[xhat2;yhat2;f0]);
        ytilde1 = dot1 / dot2 * (u(4:6)'*[xhat2;yhat2;f0]);
        
        xtilde2 = dot1 / dot2 * (u([1,4,7])'*[xhat1;yhat1;f0]);
        ytilde2 = dot1 / dot2 * (u([2,5,8])'*[xhat1;yhat1;f0]);
        
        E = (xtilde1*xtilde1+xtilde2*xtilde2+ytilde1*ytilde1+ytilde2*ytilde2)/f02;
        
        if abs(E-E0) < 1e-6
            break;
        end
        xhat1 = xhat1 - xtilde1;
        yhat1 = yhat1 - ytilde1;
        xhat2 = xhat2 - xtilde2;
        yhat2 = yhat2 - ytilde2;
        
        iter = iter + 1;
    end
    
    xu1 = [xhat1;yhat1;1];
    xu2 = [xhat2;yhat2;1];
    
    P = triangulationEigen(xu1,P1,xu2,P2);
end