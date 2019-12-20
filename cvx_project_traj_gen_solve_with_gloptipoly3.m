function cvx_project_traj_gen_solve_with_gloptipoly3
    clear all;clc;close all;
    dims = 2;
    %% test sample, Nx2
    pWpts = [[1 1]; ...
              [3 2]; ...
              [5 4]];
    vWpts = [[1 1]; ...
             [inf inf]; ...
             [0 0]];
    aWpts = [[1 0]; ...
             [inf inf]; ...
             [0 0]];
    
    vmax = 6;% 5m/s
    amax = 5;% 6.8 m/s^2
         
    %% pre-allocated timestamp, may not be reasonable
    t = 0:1:size(pWpts,1)-1;
    order = 5;

    [polyCoeffs] = trajectoryGenerator(pWpts, vWpts, aWpts, t, dims, order);
    
    visulization(pWpts, polyCoeffs, t, order, vmax, amax);
    
    lambda = 0;
    topt = [0];
    segments = size(pWpts,1)-1;
    numCoeff = order + 1;

    %% optimize t, in this case, x = [1 t t^2 t^3], do optimization for each segment
    for i = 1:segments%% each segment
        nopt = 4;
        tseg = -1e6;
    
        mpol x1
        
        %% firstly reformulate cost function
        polyx = polyCoeffs((i-1)*numCoeff+1:i*numCoeff,1);
        Cx = costPoly(polyx, x1, topt(end));
        polyy = polyCoeffs((i-1)*numCoeff+1:i*numCoeff,2);
        Cy = costPoly(polyy, x1, topt(end));
        Cxy = Cx + Cy + lambda * x1^2;

        K = [];
        poly = [polyx polyy];
        %% then do the same for position, velocity, acceleration constaints
        for j = 1:dims%% each dimension
            p0 = poly(1,j);
            p1 = poly(2,j);
            p2 = poly(3,j);
            p3 = poly(4,j);
            p4 = poly(5,j);
            p5 = poly(6,j);
            
            if pWpts(i+1,j) ~= inf
                K1 = p1*x1+p2*x1^2+p3*x1^3+p4*x1^4+p5*x1^5 <= pWpts(i+1,j)-p0;
                K2 = p1*x1+p2*x1^2+p3*x1^3+p4*x1^4+p5*x1^5 >= pWpts(i+1,j)-p0;
                K = [K, K1, K2];
            end
            if vWpts(i+1,j) ~= inf
                K1 = 2*p2*x1+3*p3*x1^2+4*p4*x1^3+5*p5*x1^4 <= vWpts(i+1,j)-p1;
                K2 = 2*p2*x1+3*p3*x1^2+4*p4*x1^3+5*p5*x1^4 >= vWpts(i+1,j)-p1;
                K = [K, K1, K2];
            end
            if aWpts(i+1,j) ~= inf
                K1 = 6*p3*x1+12*p4*x1^2+20*p5*x1^3 <= aWpts(i+1,j)-2*p2;
                K2 = 6*p3*x1+12*p4*x1^2+20*p5*x1^3 >= aWpts(i+1,j)-2*p2;
                K = [K, K1, K2];
            end
            
            %% inequality constraints
            K = [K, 2*p2*x1+3*p3*x1^2+4*p4*x1^3+5*p5*x1^4 <= vmax-p1, 2*p2*x1+3*p3*x1^2+4*p4*x1^3+5*p5*x1^4 >= -vmax-p1];
            K = [K, 6*p3*x1+12*p4*x1^2+20*p5*x1^3 <= amax-2*p2, 6*p3*x1+12*p4*x1^2+20*p5*x1^3 >= -amax-2*p2];            
        end
        Prob = msdp(min(Cxy), K);
        [status,obj] = msol(Prob);
        xsol = meas;
        %% parsing solution
        tt = double(xsol);
        if tseg < tt
            tseg = tt;
        end            
        topt = [topt;topt(end)+tt];
    end
    t
    topt
    [polyCoeffsNew] = trajectoryGenerator(pWpts, vWpts, aWpts, topt, dims, order);
    visulization(pWpts, polyCoeffsNew, topt', order, vmax, amax);
end

function [polyCoeffs] = trajectoryGenerator(pWpts, vWpts, aWpts, t, dims, order)
    %% order of trajectory 5
    numCoeff = order + 1;
    segments = size(pWpts,1)-1;
    polyCoeffs = zeros(segments*numCoeff,dims);

    %% cost matrix
    Q = zeros(numCoeff*segments, numCoeff*segments);
    
    %% equality constraint
    maxNumCons = numCoeff * segments;
    alreadySetCons = 0;
    %% traverse waypts
    for i=1:size(pWpts,1)
        if i == 1 || i == size(pWpts,1)
            alreadySetCons = alreadySetCons + double(pWpts(i,1)~=inf);
            alreadySetCons = alreadySetCons + double(vWpts(i,1)~=inf);
            alreadySetCons = alreadySetCons + double(aWpts(i,1)~=inf);            
        else
            alreadySetCons = alreadySetCons + 2*double(pWpts(i,1)~=inf);
            alreadySetCons = alreadySetCons + 2*double(vWpts(i,1)~=inf);
            alreadySetCons = alreadySetCons + 2*double(aWpts(i,1)~=inf);
        end
    end
    solveAxb = 0;
    if alreadySetCons > maxNumCons
        disp('Error, too many constaints!');
    elseif alreadySetCons == maxNumCons
        disp('Solve Ax=b!');
        solveAxb = 1;
    else
        %% continuity constraints
        for i=2:size(pWpts,1)-1
            if vWpts(i,1) == inf
                alreadySetCons = alreadySetCons + 1;
            end
            if aWpts(i,1) == inf
                alreadySetCons = alreadySetCons + 1;
            end
        end
    end
    %% cnstruct A
    Aeq = zeros(alreadySetCons,numCoeff*segments);
    beqs = zeros(alreadySetCons,dims);
    k = 1;
    for i=1:segments
        t0 = t(i);
        t1 = t(i+1);
        %% start point
        if pWpts(i,1) ~= inf
            Aeq(k,(i-1)*numCoeff+1:i*numCoeff) = [1 t0 t0^2 t0^3 t0^4 t0^5];
            beqs(k,:) = pWpts(i,:);
            k = k + 1;
        end
        if vWpts(i,1) ~= inf
            Aeq(k,(i-1)*numCoeff+1:i*numCoeff) = [0 1 2*t0 3*t0^2 4*t0^3 5*t0^4];
            beqs(k,:) = vWpts(i,:);
            k = k + 1;
        end
        if aWpts(i,1) ~= inf
            Aeq(k,(i-1)*numCoeff+1:i*numCoeff) = [0 0 2 6*t0 12*t0^2 20*t0^3];
            beqs(k,:) = aWpts(i,:);
            k = k + 1;
        end
        %% end point
        if pWpts(i+1,1) ~= inf
            Aeq(k,(i-1)*numCoeff+1:i*numCoeff) = [1 t1 t1^2 t1^3 t1^4 t1^5];
            beqs(k,:) = pWpts(i+1,:);
            k = k + 1;
        end
        if vWpts(i+1,1) ~= inf
            Aeq(k,(i-1)*numCoeff+1:i*numCoeff) = [0 1 2*t1 3*t1^2 4*t1^3 5*t1^4];
            beqs(k,:) = vWpts(i+1,:);
            k = k + 1;
        end
        if aWpts(i+1,1) ~= inf
            Aeq(k,(i-1)*numCoeff+1:i*numCoeff) = [0 0 2 6*t1 12*t1^2 20*t1^3];
            beqs(k,:) = aWpts(i+1,:);
            k = k + 1;
        end
    end
    %% continuity
    for i=2:size(pWpts,1)-1
        t0 = t(i);
        t1 = t(i);
        if vWpts(i,1) == inf
            Aeq(k,(i-2)*numCoeff+1:(i)*numCoeff) = [0 1 2*t0 3*t0^2 4*t0^3 5*t0^4 -0 -1 -2*t1 -3*t1^2 -4*t1^3 -5*t1^4];
            beqs(k,:) = zeros(1,dims);
            k = k + 1;
        end
        if aWpts(i,1) == inf
            Aeq(k,(i-2)*numCoeff+1:(i)*numCoeff) = [0 0 2 6*t0 12*t0^2 20*t0^3 -0 -0 -2 -6*t1 -12*t1^2 -20*t1^3];
            beqs(k,:) = zeros(1,dims);
            k = k + 1;
        end
    end
    
    %% Q
    for i = 1:segments
        Qs = constructQ(t(i), t(i+1));
        Q((i-1)*numCoeff+1:(i)*numCoeff, (i-1)*numCoeff+1:(i)*numCoeff) = Qs;
    end
     
    %% solve equality constraint convex optimization problem
    %% here I tried two solutions
    tic
    KKT = [Q Aeq';Aeq zeros(size(Q,1)+size(Aeq,1)-size(Aeq,2))];
    xsol = inv(KKT)*[zeros(size(Q,1),1);beqs(:,1)];
    ysol = inv(KKT)*[zeros(size(Q,1),1);beqs(:,2)];
    toc
%     tic
%     xsol = quadprog(Q,[],[],[],Aeq,beqs(:,1));
%     ysol = quadprog(Q,[],[],[],Aeq,beqs(:,2));
%     toc

    
%     tic
%     xsol = solveViaQR(Q,Aeq,beqs(:,1));
% %     toc
%     ysol = solveViaQR(Q,Aeq,beqs(:,2));
    polyCoeffs(:,1) = xsol(1:size(Q,1));
    polyCoeffs(:,2) = ysol(1:size(Q,1));
end

function p = costPoly(poly, x, ts)
    p3 = poly(4);
    p4 = poly(5);
    p5 = poly(6);

    a1 = 36*p3^2;
    a2 = 144*p3*p4;
    a3 = (192*p4^2 + 240*p3*p5);
    a4 = 720*p4*p5;
    a5 = 720*p5^2;
    
    p = a1*x+a2*x^2+a3*x^3+a4*x^4+a5*x^5 - (a1*ts+a2*ts^2+a3*ts^3+a4*ts^4+a5*ts^5);
end

function sol = solveViaQR(Q,A,b)
    [Qt Rt] = qr(A');
    p = size(A,1);
    n = size(A,2);
    Q1 = Qt(1:n,1:p);
    Q2 = Qt(1:n,p+1:end);
    R = Rt(1:p,1:p);
    x0 = A\((R'*Q1'*Q1*inv(R)')*b);
    Q11 = Q2'*Q*Q2;
    Q12 = Q*Q2;
    z = -(Q11)\(Q12'*x0);
    sol = x0 + Q2*z;
end

function Qs = constructQ(ts, te)
    syms t real
    assume(t,'positive');
    f3d = [0 0 0 6 24*t 60*t^2];
    Q1 = f3d'*f3d;
    Q2 = int(Q1,t);
    Qs = subs(Q2,te) - subs(Q2,ts);
end

function visulization(pWpts, polyCoeffs, t, order, vmax, amax)
    figure
    plot(pWpts(:,1),pWpts(:,2),'ro','MarkerSize',5);
    hold on;
    grid on;
    pts = [];
    vts = [];
    ats = [];

    numCoeff = order+1;
    for i = 1:size(t,2)-1
        ts = linspace(t(i),t(i+1),1000)';
        polyx = polyCoeffs((i-1)*numCoeff+1:i*numCoeff,1);
        polyy = polyCoeffs((i-1)*numCoeff+1:i*numCoeff,2);
        xsamples = polyx(1) + polyx(2).*ts + polyx(3).*ts.^2 + polyx(4).*ts.^3 + polyx(5).*ts.^4 + polyx(6).*ts.^5;
        ysamples = polyy(1) + polyy(2).*ts + polyy(3).*ts.^2 + polyy(4).*ts.^3 + polyy(5).*ts.^4 + polyy(6).*ts.^5;
        
        vxsamples = polyx(2) + polyx(3).*ts.*2 + polyx(4).*3.*ts.^2 + polyx(5).*4.*ts.^3 + polyx(6).*5.*ts.^4;
        vysamples = polyy(2) + polyy(3).*ts.*2 + polyy(4).*3.*ts.^2 + polyy(5).*4.*ts.^3 + polyy(6).*5.*ts.^4;
        
        axsamples = polyx(3).*2 + polyx(4).*6.*ts.^1 + polyx(5).*12.*ts.^2 + polyx(6).*20.*ts.^3;
        aysamples = polyy(3).*2 + polyy(4).*6.*ts.^1 + polyy(5).*12.*ts.^2 + polyy(6).*20.*ts.^3;
        
        pts = [pts;[xsamples ysamples]];
        vts = [vts;[vxsamples vysamples]];
        ats = [ats;[axsamples aysamples]];
    end
    
    plot(pts(:,1),pts(:,2),'b-');
    title('Trajectory');
    figure
    subplot(2,1,1);
    plot(vts(:,1),'b-');grid on;hold on;
    plot(ones(numel(vts(:,1)),1).*vmax,'r.-');
    plot(ones(numel(vts(:,1)),1).*-vmax,'r.-');
    legend('Generated','Maximum');
    title('Velocity');
    subplot(2,1,2);
    plot(vts(:,2),'b-');grid on;hold on;
    plot(ones(numel(vts(:,2)),1).*vmax,'r.-');
    plot(ones(numel(vts(:,2)),1).*-vmax,'r.-');
    legend('Generated','Maximum');
    figure    
    subplot(2,1,1);
    plot(ats(:,1),'b-');grid on;hold on;
    plot(ones(numel(ats(:,1)),1).*amax,'r.-');
    plot(ones(numel(ats(:,1)),1).*-amax,'r.-');
    legend('Generated','Maximum');
    title('Acceleration');
    subplot(2,1,2);
    plot(ats(:,2),'b-');grid on;hold on;
    plot(ones(numel(ats(:,2)),1).*amax,'r.-');
    plot(ones(numel(ats(:,2)),1).*-amax,'r.-');
    legend('Generated','Maximum');
end


