function cvx_project_traj_gen_mix
    clear all;clc;close all;
    dims = 2;
    %% test sample, Nx2
    pWpts = [[1 1]; ...
              [3 2]; ...
              [5 2]];
    vWpts = [[1 1]; ...
             [inf inf]; ...
             [0 0]];
    aWpts = [[1 0]; ...
             [inf inf]; ...
             [0 0]];
    
    vmax = 6;% 5m/s
    amax = 4;% 6.8 m/s^2
         
    %% pre-allocated timestamp, may not be reasonable
    t = 0:1:size(pWpts,1)-1;
    order = 5;

    [polyCoeffs] = trajectoryGenerator(pWpts, vWpts, aWpts, t, dims, order);
    
    visulization(pWpts, polyCoeffs, t, order, vmax, amax);
    
    lambda = 100;
    topt = [0];
    segments = size(pWpts,1)-1;
    numCoeff = order + 1;
    margin = 0.0;

    %% optimize t, in this case, x = [1 t t^2 t^3], do optimization for each segment
    for i = 1:segments%% each segment
        nopt = 4;
        tseg = -1e6;
    
        %% firstly reformulate cost function
        Cx = zeros(nopt,nopt);
        polyx = polyCoeffs((i-1)*numCoeff+1:i*numCoeff,1);
        Cx = costC(polyx);

        Cy = zeros(nopt,nopt);
        polyy = polyCoeffs((i-1)*numCoeff+1:i*numCoeff,2);
        Cy = costC(polyy);
        Cxy = Cx+Cy+[0 0 0 0;1 0 0 0;1 0 0 0;1 0 0 0].*lambda;

        Aopt = [];
        bopt = [];
        k = 1;

        poly = [polyx polyy];

        %% then do the same for position, velocity, acceleration constaints
        for j = 1:dims%% each dimension
            if pWpts(i,j) ~= inf
                C1 = makeC(poly(:,j), 'p',0);
                Aopt(:,:,k) = C1;
                bopt(k) = pWpts(i,j)-margin;
                Aopt(:,:,k+1) = -C1;
                bopt(k+1) = -pWpts(i,j)+margin;
                k = k+2;
            end
            if vWpts(i,j) ~= inf
                C1 = makeC(poly(:,j), 'v',0);
                Aopt(:,:,k) = C1;
                bopt(k) = vWpts(i,j)-margin;
                Aopt(:,:,k+1) = -C1;
                bopt(k+1) = -vWpts(i,j)-margin;
                k = k+2;
            end
            if aWpts(i,j) ~= inf
                C1 = makeC(poly(:,j), 'a',0);
                Aopt(:,:,k) = C1;
                bopt(k) = aWpts(i,j)-margin;
                Aopt(:,:,k+1) = -C1;
                bopt(k+1) = -aWpts(i,j)-margin;
                k = k+2;
            end
            if pWpts(i+1,j) ~= inf
                C1 = makeC(poly(:,j), 'p',1);
                Aopt(:,:,k) = C1;
                bopt(k) = pWpts(i+1,j)-margin;
                Aopt(:,:,k+1) = -C1;
                bopt(k+1) = -pWpts(i+1,j)-margin;
                k = k+2;
            end
            if vWpts(i+1,j) ~= inf
                C1 = makeC(poly(:,j), 'v',1);
                Aopt(:,:,k) = C1;
                bopt(k) = vWpts(i+1,j)-margin;
                Aopt(:,:,k+1) = -C1;
                bopt(k+1) = -vWpts(i+1,j)-margin;
                k = k+2;
            end
            if aWpts(i+1,j) ~= inf
                C1 = makeC(poly(:,j), 'a',1);
                Aopt(:,:,k) = C1;
                bopt(k) = aWpts(i+1,j)-margin;
                Aopt(:,:,k+1) = -C1;
                bopt(k+1) = -aWpts(i+1,j)-margin;
                k = k+2;
            end
            numeq = k;
            %% inequality constraints
            Cvs = makeC(poly(:,j), 'v',0);
            Aopt(:,:,k) = -Cvs;bopt(k) = vmax;
            Aopt(:,:,k+1) = Cvs;bopt(k+1) = vmax;
            k = k + 2;
            Cve = makeC(poly(:,j), 'v',1);
            Aopt(:,:,k) = -Cve;bopt(k) = vmax;
            Aopt(:,:,k+1) = Cve;bopt(k+1) = vmax;
            k = k + 2;
            Cas = makeC(poly(:,j), 'a',0);
            Aopt(:,:,k) = -Cas;bopt(k) = amax;
            Aopt(:,:,k+1) = Cas;bopt(k+1) = amax;
            k = k + 2;
            Cae = makeC(poly(:,j), 'a',1);
            Aopt(:,:,k) = -Cae;bopt(k) = amax;
            Aopt(:,:,k+1) = Cae;bopt(k+1) = amax;
            k = k + 2;
        end
            
        %% cvx_routine
        cvx_begin
            variable X(nopt, nopt) nonnegative hankel
            minimize( trace(Cxy*X) );
            subject to 
%                 for ii = 1:numeq-1
%                     trace(Aopt(:,:,ii)*X) >= bopt(ii);
%                 end
                for ii = 1:k-1
                    trace(Aopt(:,:,ii)*X) <= bopt(ii);
                end
                X == semidefinite(nopt)
                X(1,1) == 1
                X(2,1) >= 0
                X(3,1) >= 0
                X(4,1) >= 0
                X(4,2) >= 0
                X(4,3) >= 0
                X(4,4) >= 0
        cvx_end

        %% here, I just take the maximum eigvalue corresponding solution
        [eigvec, eigval] = eig(X);
        maxEig = eigval(end,end);
        correspondingEvec = eigvec(:,end);
        xsol = sqrt(maxEig).*correspondingEvec;
        rank(X)
        %% parsing solution
        tt = abs(xsol(4))^(1/3);
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
        t0 = 0;
        t1 = t(i+1)-t(i);
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
        t0 = t(i)-t(i-1);
        t1 = 0;
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
        ts = t(i+1)-t(i);
        Qs = constructQ(ts);
        Q((i-1)*numCoeff+1:(i)*numCoeff, (i-1)*numCoeff+1:(i)*numCoeff) = Qs;
    end
    
    %% constuct C
%     A1 = zeros(numCoeff*segments,numCoeff*segments);
%     for i=1:segments
%         t0 = 0;
%         t1 = t(i+1)-t(i);
%         %% start point
%         A1((i-1)*6+1,(i-1)*numCoeff+1:i*numCoeff) = [1 t0 t0^2 t0^3 t0^4 t0^5];
%         A1((i-1)*6+2,(i-1)*numCoeff+1:i*numCoeff) = [0 1 2*t0 3*t0^2 4*t0^3 5*t0^4];
%         A1((i-1)*6+3,(i-1)*numCoeff+1:i*numCoeff) = [0 0 2 6*t0 12*t0^2 20*t0^3];
%         %% end point
%         A1((i-1)*6+4,(i-1)*numCoeff+1:i*numCoeff) = [1 t1 t1^2 t1^3 t1^4 t1^5];
%         A1((i-1)*6+5,(i-1)*numCoeff+1:i*numCoeff) = [0 1 2*t1 3*t1^2 4*t1^3 5*t1^4];
%         A1((i-1)*6+6,(i-1)*numCoeff+1:i*numCoeff) = [0 0 2 6*t1 12*t1^2 20*t1^3];
%     end
%     C = [1 0 0 0 0 0 0 0 0; ...
%          0 1 0 0 0 0 0 0 0; ...
%          0 0 1 0 0 0 0 0 0; ...
%          0 0 0 1 0 0 0 0 0; ...
%          0 0 0 0 0 0 0 1 0; ...
%          0 0 0 0 0 0 0 0 1; ...
%          0 0 0 1 0 0 0 0 0; ...
%          0 0 0 0 0 0 0 1 0; ...
%          0 0 0 0 0 0 0 0 1; ...
%          0 0 0 0 1 0 0 0 0; ...
%          0 0 0 0 0 1 0 0 0; ...
%          0 0 0 0 0 0 1 0 0];
%      tic
%      R = C'*inv(A1)'*Q*inv(A1)*C;
%      Rpp = R(8:9,8:9);
%      Rfp = R(1:7,8:9);
%      dp = -inv(Rpp)*Rfp'*[1 1 1 3 5 0 0]';
%      inv(A1)*C*[[1 1 1 3 5 0 0]';dp];
%      toc
     
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

function C = makeC(poly, mode, id)
    p0 = poly(1);
    p1 = poly(2);
    p2 = poly(3);
    p3 = poly(4);
    p4 = poly(5);
    p5 = poly(6);
    if mode == 'p'
        if id == 0
            C = [p0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0];
        elseif id == 1
            C = [p0 0 0 0; ...
                 p1 0 0 0; ...
                 p2 0 0 0; ...
                 p3 p4 p5 0];
        end
    elseif mode == 'v'
        if id == 0
            C = [p1 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0];
        elseif id == 1
            C = [p1 0 0 0; ...
                 2*p2 0 0 0; ...
                 3*p3 0 0 0; ...
                 4*p4 5*p5 0 0];
        end
    elseif mode == 'a'
        if id == 0
            C = [2*p2 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0];
        elseif id == 1
            C = [2*p2 0 0 0; ...
                 6*p3 0 0 0; ...
                 12*p4 0 0 0; ...
                 20*p5 0 0 0];
        end
    end
end

function C = costC(poly)
    p3 = poly(4);
    p4 = poly(5);
    p5 = poly(6);

    a1 = 36*p3^2;
    a2 = 144*p3*p4;
    a3 = (192*p4^2 + 240*p3*p5);
    a4 = 720*p4*p5;
    a5 = 720*p5^2;
    
    C = [0 0 0 0;...
         a1 0 0 0;...
         a2 0 0 0;...
         a3 a4 a5 0];
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

function feval = evaluateFunc()
end

function Qs = constructQ(ts)
    syms t real
    assume(t,'positive');
    f3d = [0 0 0 6 24*t 60*t^2];
    Q1 = f3d'*f3d;
    Q2 = int(Q1,t);
    Qs = subs(Q2,ts);
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
        ts = linspace(0,t(i+1)-t(i),1000)';
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


