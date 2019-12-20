function [polyCoeffs,realt]=cvx_project_traj_gen_solver(p,initv,inita,endv,enda,maxv,maxa)
    dims = 2;
    N = size(p,1);
    pWpts = zeros(N,2);
    vWpts = zeros(N,2);
    aWpts = zeros(N,2);
     
    pWpts = p;
    vWpts(1,1) = initv(1);
    vWpts(1,2) = initv(2);
    vWpts(end,1) = endv(1);
    vWpts(end,2) = endv(2);
    vWpts(2:end-1,:) = inf;

    aWpts(1,1) = inita(1);
    aWpts(1,2) = inita(2);
    aWpts(end,1) = enda(1);
    aWpts(end,2) = enda(2);
    aWpts(2:end-1,:) = inf;

    %% test sample, Nx2
%     pWpts = [[1 1]; ...
%               [3 2]; ...
%               [5 4]];
%     vWpts = [[1 1]; ...
%              [inf inf]; ...
%              [0 0]];
%     aWpts = [[1 0]; ...
%              [inf inf]; ...
%              [0 0]];
    
    vmax = maxv;% 5m/s
    amax = maxa * 0.99;% 6.8 m/s^2
    
    t = linspace(0,size(pWpts,1)-1,size(pWpts,1));
    %t = [0 2 4];
    order = 5;

    %% 1. pre-allocated timestamp, may not be reasonable
    [polyCoeffs] = trajectoryGenerator(pWpts, vWpts, aWpts, t, dims, order);
%     visulization(pWpts, polyCoeffs, t, order, vmax, amax);    
    
    origS = [1;1];
    realt = t;

    while 1
        segments = size(pWpts,1)-1;
        numCoeff = order + 1;
        replan = 0;
        %% 2. find extreme point
        vxext = {};vyext = {};
        axext = {};ayext = {};
        for i=1:segments
            scalar = 1./(realt(i+1)-realt(i));
            
            polyx = polyCoeffs((i-1)*numCoeff+1:i*numCoeff,1);
            polyy = polyCoeffs((i-1)*numCoeff+1:i*numCoeff,2);

            vxpoly = scalar.*[5*polyx(6) 4*polyx(5) 3*polyx(4) 2*polyx(3) polyx(2)];
            vypoly = scalar.*[5*polyy(6) 4*polyy(5) 3*polyy(4) 2*polyy(3) polyy(2)];
            axpoly = scalar^2.*[20*polyx(6) 12*polyx(5) 6*polyx(4) 2*polyx(3)];
            aypoly = scalar^2.*[20*polyy(6) 12*polyy(5) 6*polyy(4) 2*polyy(3)];

            vdxpoly = scalar^2.*[20*polyx(6) 12*polyx(5) 6*polyx(4) 2*polyx(3)];
            vdypoly = scalar^2.*[20*polyy(6) 12*polyy(5) 6*polyy(4) 2*polyy(3)];

            adxpoly = scalar^3.*[60*polyx(6) 24*polyx(5) 6*polyx(4)];
            adypoly = scalar^3.*[60*polyy(6) 24*polyy(5) 6*polyy(4)];

            xroots = roots(vdxpoly);
            yroots = roots(vdypoly);
            idx = isreal(xroots) & real(xroots) > 0 & real(xroots) < 1;
            xroots = xroots(idx);
            idy = isreal(yroots) & real(yroots) > 0 & real(yroots) < 1;
            yroots = yroots(idy);
            xroots = [xroots;0;1];
            yroots = [yroots;0;1];

            if ~isempty(xroots)
                vxmaxmin = polyval(vxpoly,xroots);
            end
            if ~isempty(yroots)
                vymaxmin = polyval(vypoly,yroots);
            end

            xroots = roots(adxpoly);
            yroots = roots(adypoly);
            idx = isreal(xroots) & real(xroots) > 0 & real(xroots) < 1;
            xroots = xroots(idx);
            idy = isreal(yroots) & real(yroots) > 0 & real(yroots) < 1;
            yroots = yroots(idy);
            xroots = [xroots;0;1];
            yroots = [yroots;0;1];
            if ~isempty(xroots)
                axmaxmin = polyval(axpoly,xroots);
            end
            if ~isempty(yroots)
                aymaxmin = polyval(aypoly,yroots);
            end
            vxext{i} = vxmaxmin;
            vyext{i} = vymaxmin;
            axext{i} = axmaxmin;
            ayext{i} = aymaxmin;
            
            r1 = isempty(find(abs(vxmaxmin) > vmax)) && ...
                     isempty(find(abs(vymaxmin) > vmax)) && ...
                     isempty(find(abs(axmaxmin) > amax)) && ...
                     isempty(find(abs(aymaxmin) > amax));
            replan = replan + r1;
        end

        lambda = 0;
        
        if replan == segments
            break;
        end
        
        %% try to find an linear mapping scalar for each segment, s1, s2, ..., sn
        %% the scalar is used to map from real world timestamp to virtual timestamp
        scalars = [];
        
        useCVX = 1;
        if useCVX == 0
            mpol s 2
            %% firstly reformulate cost function
            g0 = s(1)^2-2*1*s(1)+s(2)^2-2*1*s(2)+2;
            %% then do the same for position, velocity, acceleration constaints
            K = [];
            %% vx
            for j = 1:segments
                for k = 1:size(vxext{j},1)
                    K = [K, vxext{j}(k)*s(j)<=vmax, vxext{j}(k)*s(j)>=-vmax];
                end
            end
            %% vy
            for j = 1:segments
                for k = 1:size(vyext{j},1)
                    K = [K, vyext{j}(k)*s(j)<=vmax, vyext{j}(k)*s(j)>=-vmax];
                end
            end
            %% ax
            for j = 1:segments
                for k = 1:size(axext{j},1)
                    K = [K, axext{j}(k)*s(j)^2<=amax, axext{j}(k)*s(j)^2>=-amax];
                end
            end
            %% ay
            for j = 1:segments
                for k = 1:size(ayext{j},1)
                    K = [K, ayext{j}(k)*s(j)^2<=amax, ayext{j}(k)*s(j)^2>=-amax];
                end
            end
            K = [K, s(1)>=0, s(2)>=0];
            Prob = msdp(min(g0), K);
            [status,obj] = msol(Prob);
            xsol = meas;
            %% parsing solution
            tt = double(xsol);
            scalars = [scalars;tt];
        else
            %% USE CVX, model this as a SDP problem with SDR
            C0 = [1 -2;0 1];
            Aopt = [];
            bopt = [];
            nopt = 2;
            %% vx
            for j = 1:segments
                kk = 1;
                for k = 1:size(vxext{j},1)
                    Aopt(:,:,kk) = [0 vxext{j}(k); 0 0];
                    bopt(kk) = vmax* 0.8;
                    Aopt(:,:,kk+1) = [0 -vxext{j}(k); 0 0];
                    bopt(kk+1) = vmax* 0.8;
                    kk = kk + 2;
                end
                for k = 1:size(vyext{j},1)
                    Aopt(:,:,kk) = [0 vyext{j}(k); 0 0];
                    bopt(kk) = vmax* 0.8;
                    Aopt(:,:,kk+1) = [0 -vyext{j}(k); 0 0];
                    bopt(kk+1) = vmax* 0.8;
                    kk = kk + 2;
                end
                for k = 1:size(axext{j},1)
                    Aopt(:,:,kk) = [0 0; 0 axext{j}(k)];
                    bopt(kk) = amax * 0.8;
                    Aopt(:,:,kk+1) = [0 0; 0 -axext{j}(k)];
                    bopt(kk+1) = amax * 0.8;
                    kk = kk + 2;
                end
                for k = 1:size(ayext{j},1)
                    Aopt(:,:,kk) = [0 0; 0 ayext{j}(k)];
                    bopt(kk) = amax * 0.8;
                    Aopt(:,:,kk+1) = [0 0; 0 -ayext{j}(k)];
                    bopt(kk+1) = amax * 0.8;
                    kk = kk + 2;
                end
                %% cvx_routine
                cvx_begin quiet
                    variable X(nopt, nopt) nonnegative hankel
                    minimize( trace(C0*X) );
                    subject to 
                        for ii = 1:kk-1
                            trace(Aopt(:,:,ii)*X) <= bopt(ii);
                        end
                        X == semidefinite(nopt)
                        X(1,1) == 1
                        X(1,2) >= 0
                        X(2,2) >= 0
                cvx_end
                %% here, I just take the maximum eigvalue corresponding solution
                [eigvec, eigval] = eig(X);
                maxEig = eigval(end,end);
                correspondingEvec = eigvec(:,end);
                xsol = sqrt(maxEig).*correspondingEvec;
                rank(X)
                %% parsing solution
                tt = double(xsol(2)/xsol(1));
                scalars = [scalars;tt];
            end
        end
        ttmp = realt;
        for i = 2:numel(t)
            realt(i) = realt(i-1) + (ttmp(i)-ttmp(i-1)) / scalars(i-1);
        end
        
        origS = scalars;
        
        [polyCoeffs] = trajectoryGenerator(pWpts, vWpts, aWpts, realt, dims, order);
    end
%     visulization(pWpts, polyCoeffs, realt, order, vmax, amax);
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
        t1 = 1;%t(i+1)-t(i);
        scalar = 1./(t(i+1)-t(i));
        %% start point
        if pWpts(i,1) ~= inf
            Aeq(k,(i-1)*numCoeff+1:i*numCoeff) = [1 t0 t0^2 t0^3 t0^4 t0^5];
            beqs(k,:) = pWpts(i,:);
            k = k + 1;
        end
        if vWpts(i,1) ~= inf
            Aeq(k,(i-1)*numCoeff+1:i*numCoeff) = scalar.*[0 1 2*t0 3*t0^2 4*t0^3 5*t0^4];
            beqs(k,:) = vWpts(i,:);
            k = k + 1;
        end
        if aWpts(i,1) ~= inf
            Aeq(k,(i-1)*numCoeff+1:i*numCoeff) = scalar^2.*[0 0 2 6*t0 12*t0^2 20*t0^3];
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
            Aeq(k,(i-1)*numCoeff+1:i*numCoeff) = scalar.*[0 1 2*t1 3*t1^2 4*t1^3 5*t1^4];
            beqs(k,:) = vWpts(i+1,:);
            k = k + 1;
        end
        if aWpts(i+1,1) ~= inf
            Aeq(k,(i-1)*numCoeff+1:i*numCoeff) = scalar^2.*[0 0 2 6*t1 12*t1^2 20*t1^3];
            beqs(k,:) = aWpts(i+1,:);
            k = k + 1;
        end
    end
    %% continuity
    for i=2:size(pWpts,1)-1
        t0 = 1;
        t1 = 0;
        s1 = 1./(t(i)-t(i-1));
        s2 = 1./(t(i+1)-t(i));
        if vWpts(i,1) == inf
            Aeq(k,(i-2)*numCoeff+1:(i-1)*numCoeff) = s1.*[0 1 2*t0 3*t0^2 4*t0^3 5*t0^4];
            Aeq(k,(i-1)*numCoeff+1:(i)*numCoeff) = s2.*[-0 -1 -2*t1 -3*t1^2 -4*t1^3 -5*t1^4];
            beqs(k,:) = zeros(1,dims);
            k = k + 1;
        end
        if aWpts(i,1) == inf
            Aeq(k,(i-2)*numCoeff+1:(i-1)*numCoeff) = s1^2.*[0 0 2 6*t0 12*t0^2 20*t0^3];
            Aeq(k,(i-1)*numCoeff+1:(i)*numCoeff) = s2^2.*[-0 -0 -2 -6*t1 -12*t1^2 -20*t1^3];
            beqs(k,:) = zeros(1,dims);
            k = k + 1;
        end
    end
    
    %% Q
    for i = 1:segments
        Qs = constructQ(0, t(i+1)-t(i));
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
    scalar = 1/te;
    f3d = scalar^3.*[0 0 0 6 24*t 60*t^2];
    Q1 = f3d'*f3d;
    Q2 = int(Q1,t);
    Qs = subs(Q2,1) - subs(Q2,ts);
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
    tss = [];
    for i = 1:size(t,2)-1
        scalar = 1./(t(i+1)-t(i));
        ts = linspace(0,1,1000)';
        polyx = polyCoeffs((i-1)*numCoeff+1:i*numCoeff,1);
        polyy = polyCoeffs((i-1)*numCoeff+1:i*numCoeff,2);
        xsamples = polyx(1) + polyx(2).*ts + polyx(3).*ts.^2 + polyx(4).*ts.^3 + polyx(5).*ts.^4 + polyx(6).*ts.^5;
        ysamples = polyy(1) + polyy(2).*ts + polyy(3).*ts.^2 + polyy(4).*ts.^3 + polyy(5).*ts.^4 + polyy(6).*ts.^5;
        
        vxsamples = polyx(2) + polyx(3).*ts.*2 + polyx(4).*3.*ts.^2 + polyx(5).*4.*ts.^3 + polyx(6).*5.*ts.^4;
        vysamples = polyy(2) + polyy(3).*ts.*2 + polyy(4).*3.*ts.^2 + polyy(5).*4.*ts.^3 + polyy(6).*5.*ts.^4;
        
        axsamples = polyx(3).*2 + polyx(4).*6.*ts.^1 + polyx(5).*12.*ts.^2 + polyx(6).*20.*ts.^3;
        aysamples = polyy(3).*2 + polyy(4).*6.*ts.^1 + polyy(5).*12.*ts.^2 + polyy(6).*20.*ts.^3;
        
        vxsamples = vxsamples.*scalar;
        vysamples = vysamples.*scalar;
        axsamples = axsamples.*(scalar^2);
        aysamples = aysamples.*(scalar^2);
        
        pts = [pts;[xsamples ysamples]];
        vts = [vts;[vxsamples vysamples]];
        ats = [ats;[axsamples aysamples]];
        
        tss = [tss;ts+(i-1)];
    end
    
    plot(pts(:,1),pts(:,2),'b-');
    title('Trajectory');
    figure
    subplot(2,1,1);
    plot(tss,vts(:,1),'b-');grid on;hold on;
    plot(tss,ones(numel(vts(:,1)),1).*vmax,'r.-');
    plot(tss,ones(numel(vts(:,1)),1).*-vmax,'r.-');
    legend('Generated','Maximum');
    title('Velocity');
    subplot(2,1,2);
    plot(tss,vts(:,2),'b-');grid on;hold on;
    plot(tss,ones(numel(vts(:,2)),1).*vmax,'r.-');
    plot(tss,ones(numel(vts(:,2)),1).*-vmax,'r.-');
    legend('Generated','Maximum');
    figure    
    subplot(2,1,1);
    plot(tss,ats(:,1),'b-');grid on;hold on;
    plot(tss,ones(numel(ats(:,1)),1).*amax,'r.-');
    plot(tss,ones(numel(ats(:,1)),1).*-amax,'r.-');
    legend('Generated','Maximum');
    title('Acceleration');
    subplot(2,1,2);
    plot(tss,ats(:,2),'b-');grid on;hold on;
    plot(tss,ones(numel(ats(:,2)),1).*amax,'r.-');
    plot(tss,ones(numel(ats(:,2)),1).*-amax,'r.-');
    legend('Generated','Maximum');
end


