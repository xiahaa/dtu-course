function [R,t] = hand_eye_calib_sol2()

    %% simulated data
    N = 10;
    T1 = cell(N,1);
    T2 = cell(N,1);
    
    dR = eul2rotm([rand(1,3).*pi./180]);
    dt = rand(3,1) * 2 - 1;
    dT = [dR dt;[0 0 0 1]];
    for i = 1:N
        R1 = eul2rotm([rand(1,3).*pi./180]);
%         R2 = dR*R1;
        t1 = rand(3,1) * 10 - 5;
%         t2 = (dR*(t1) + dt);
        T1{i} = [R1 t1;[0 0 0 1]];
        T = dT;
        T2{i} = dT*T1{i}*[T(1:3,1:3)' -T(1:3,1:3)'*T(1:3,4);[0 0 0 1]];
%         
%         test = [T(1:3,1:3)' -T(1:3,1:3)'*T(1:3,4);[0 0 0 1]]*([dR dt;[0 0 0 1]]*T1{i})
    end
    
    %%
    A = zeros(N*3,3);
    b = zeros(N*3,1);
    for i = 1:N
        Ta = T1{i};
        Tb = T2{i};
        Ra = Ta(1:3,1:3);
        Rb = Tb(1:3,1:3);
        [pra] = R2pr(Ra);
        [prb] = R2pr(Rb);
        A((i-1)*3+1:i*3,1:3) = skewm(pra+prb);
        b((i-1)*3+1:i*3) = pra - prb;
    end
    y = A\b;
    p12 = 2.*y./sqrt(1+y'*y);
    R12 = pr2R(p12);
    
    A = zeros(N*3,3);
    b = zeros(N*3,1);
    for i = 1:N
        Ta = T1{i};
        Tb = T2{i};
        Rb = Tb(1:3,1:3);
        A((i-1)*3+1:i*3,1:3) = Rb - eye(3);
        b((i-1)*3+1:i*3) = R12*Ta(1:3,4)-Tb(1:3,4);
    end
    t12 = A\b;
    
    disp(dt);
    disp(dR);
    disp(t12);
    disp(R12);
    
    [R,t] = hand_eye_cvx(T1, T2, N);
    disp(R)
    disp(t)
end

function R = pr2R(pr)
    pr2 = pr'*pr;
    R = (1-(pr2)*0.5)*eye(3)+0.5.*(pr*pr'+sqrt(4-pr2).*skewm(pr));
end

function [pr] = R2pr(R)
    format long
    theta = acos((trace(R)-1)*0.5);
    [U,S,V] = svd(R-eye(3));
    pr = V(:,end).*2.*sin(theta*0.5);
    if trace(pr2R(pr)'*R-eye(3)) > 1e-3
        pr = -V(:,end).*2.*sin(theta*0.5);
    end
    
%     disp(pr2R(pr)'*R-eye(3));
end

function s = skewm(w)
    s = [0 -w(3) w(2);w(3) 0 -w(1);-w(2) w(1) 0];
end

function [R,t] = hand_eye_cvx(T1, T2, N)
    Cs = zeros(12*N,12);
    ds = zeros(12*N,1);
    for i = 1:N
        Ta = T1{i};
        Tb = T2{i};
        Ra = Ta(1:3,1:3);
        Rb = Tb(1:3,1:3);
        ta = Ta(1:3,4);
        tb = Tb(1:3,4);
        id = (i-1)*12;
        Cs(id+1:id+9,1:9) = eye(9) - kron(Ra,Rb);
        Cs(id+10:id+12,1:9) = kron(eye(3),tb');
        Cs(id+10:id+12,10:12) = eye(3)-Ra;
        ds(id+10:id+12) = ta;
    end
    C = zeros(12,12);
    d = zeros(12,1);
    %% this should be enhanced:
    %% 1. adding constraints on rotation matrix by adding regularization
    %% 2. from socp to QP
    
    cvx_begin quiet
        variable t 
        variable x(12)
        minimize( t )
        subject to
            for i = 1:N
                id = (i-1)*12;
                C = Cs(id+1:id+12,1:12);
                d = ds(id+1:id+12);
                t >= norm( C*x-d, 2);
            end
    cvx_end
    R = [x(1) x(4) x(7); ...
         x(2) x(5) x(8); ...
         x(3) x(6) x(9)];
    R = inv(sqrtm(R*R'))*R;
    t = x(10:12);
end

function [R,t] = hand_eye_dual_quaternion(T1, T2, N)

    %% in my opinion
    % 1. R2q
    % 2. qdual;
    % 3. check scalar
    % 4. T
    % 5. svd
    % 6. compute s
    % 7. compute l1 l2
    % 8. compute q qdual
    % 9. R t.

end

function q = vec2q(v)
    q = [0;v];
end

function qdual = dualq(q, t)
    qdual = 0.5.*qprod(t,q);
end

function qp = qprod(q1,q2)
    qp(1) = q1(1)*q2(1)-q1(2:4)'*q2(2:4);
    qp(2:4) = q1(1).*q2(2:4)+q2(1).*q1(2:4)+cross(q1(2:4),q2(2:4));
end

function [d, l, theta, m] = sol_screw(R,t)
    theta = acos((trace(R)-1)*0.5);
    [U,S,V]= svd(R-eye(3));
    l = V(:,end);
    Rc = eye(3) + sin(theta).*skewm(l) + (1-cos(theta)).*(skewm(l)*skewm(l));
    if norm(diag(Rc-R)) > 1e-3
        l = -l;
    end
    d = t'*l;
    m = 0.5.*(cross(t,l)+cross(l,cross(t,l)).*cot(theta*0.5));
end
