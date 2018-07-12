function velocity_profile

pWpts = [0 0;2 1;4 5];
vWpts = [0 0;inf inf;0 0];
aWpts = [0 0;inf inf;0 0];

t = 0:0.5:1;
dims = 2;
order = 5;
vmax = 8;
amax = 6;

[polyCoeffs] = trajectoryGenerator(pWpts, vWpts, aWpts, t, dims, order, [], [], []);

visulization(pWpts, polyCoeffs, t, order, vmax, amax);


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
        
        tss = [tss;linspace(t(i),t(i+1),1000)'];
    end
    
    plot(pts(:,1),pts(:,2),'b-');
    xlabel('x: (m)','Interpreter','latex');
    ylabel('y: (m)','Interpreter','latex');
    title('Trajectory','Interpreter','latex');
    legend('Waypoints','Trajectory','Interpreter','latex');
    figure
    vx1 = subplot(2,1,1);
%     plot(tss,vts(:,1),'b-');grid on;hold on;
    c = sqrt(vts(:,1).^2);
    map = colormap(jet);
    h1 = ccplot(tss,vts(:,1),c,map);grid on;hold on;
    h2 = plot(tss,ones(numel(vts(:,1)),1).*vmax,'r--','LineWidth', 3);
    plot(tss,ones(numel(vts(:,1)),1).*-vmax,'r--','LineWidth', 3);
    legend([h1(1) h2],'Generated','Maximum','Interpreter','latex');
    xlabel('Time: (s)','Interpreter','latex');
    ylabel('$V_x: (m/s)$','Interpreter','latex');
    title('Velocity Profile','Interpreter','latex');
    subplot(2,1,2);
    c = sqrt(vts(:,2).^2);
    map = colormap(jet);
    h3 = ccplot(tss,vts(:,2),c,map);grid on;hold on;
    h4 = plot(tss,ones(numel(vts(:,2)),1).*vmax,'r--','LineWidth', 3);
    plot(tss,ones(numel(vts(:,2)),1).*-vmax,'r--','LineWidth', 3);
    patch([0.45 0.85 0.85 0.45], [max(ylim)*2 max(ylim)*2 min(ylim)*2 min(ylim)*2], [0.1 0.1 0.1], 'FaceAlpha',0.2,'EdgeColor','none');
    legend([h3(1) h4],'Generated','Maximum','Interpreter','latex');
    xlabel('Time: (s)','Interpreter','latex');
    ylabel('$V_y: (m/s)$','Interpreter','latex');
    hp1 = get(subplot(2,1,1),'Position');
    hp2 = get(subplot(2,1,2),'Position');
    cv = colorbar('Position', [hp1(1)+hp1(3)+0.01 hp2(2) 0.01 hp2(4)*2.4]);
    cv.Label.String = '$V=\sqrt{(V_x^2+V_y^2)} (m/s)$';
    cv.Label.Interpreter = 'latex';
    cv.FontSize = 12;
    caxis([0,max(sqrt(vts(:,2).^2+vts(:,1).^2))]);
    figure    
    vx2 = subplot(2,1,1);
    ca = sqrt(ats(:,1).^2);
    map1 = colormap(jet);
%     plot(tss,ats(:,1),'b-');grid on;hold on;
    h5 = ccplot(tss,ats(:,1),ca,map1);grid on;hold on;
    plot(tss,ones(numel(ats(:,1)),1).*amax,'r--', 'LineWidth', 3);
    h6 = plot(tss,ones(numel(ats(:,1)),1).*-amax,'r--','LineWidth', 3);
    
    patch([0.02 0.45 0.45 0.02], [max(ylim) max(ylim) min(ylim) min(ylim)], [0.1 0.1 0.1], 'FaceAlpha',0.2,'EdgeColor','none');
    patch([0.55 0.97 0.97 0.55], [max(ylim) max(ylim) min(ylim) min(ylim)], [0.1 0.1 0.1], 'FaceAlpha',0.2,'EdgeColor','none');

    
    legend([h5(1) h6],'Generated','Maximum','Interpreter','latex');
    xlabel('Time: (s)','Interpreter','latex');
    ylabel('$A_x: (m/s)$','Interpreter','latex');
    title('Acceleration Profile','Interpreter','latex');
    subplot(2,1,2);
    ca = sqrt(ats(:,2).^2);
%     plot(tss,ats(:,2),'b-');grid on;hold on;
    h7 = ccplot(tss,ats(:,2),ca,map1);grid on;hold on;
    h8 = plot(tss,ones(numel(ats(:,2)),1).*amax,'r--','LineWidth', 3);
    plot(tss,ones(numel(ats(:,2)),1).*-amax,'r--','LineWidth', 3);
    
    patch([0.19 0.62 0.62 0.19], [max(ylim) max(ylim) min(ylim) min(ylim)], [0.1 0.1 0.1], 'FaceAlpha',0.2,'EdgeColor','none');
    patch([0.66 0.98 0.98 0.66], [max(ylim) max(ylim) min(ylim) min(ylim)], [0.1 0.1 0.1], 'FaceAlpha',0.2,'EdgeColor','none');
    
    legend([h7(1) h8],'Generated','Maximum','Interpreter','latex');
    xlabel('Time: (s)','Interpreter','latex');
    ylabel('$A_y: (m/s)$','Interpreter','latex');
    hp3 = get(subplot(2,1,1),'Position');
    hp4 = get(subplot(2,1,2),'Position');
    cv1 = colorbar('Position', [hp3(1)+hp3(3)+0.01 hp4(2) 0.01 hp4(4)*2.4]);
    cv1.Label.String = '$A=\sqrt{(A_x^2+A_y^2)} (m/s^2)$';
    cv1.Label.Interpreter = 'latex';
    cv1.FontSize = 12;
    caxis([0,max(sqrt(ats(:,1).^2 + ats(:,2).^2))]);
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

function [polyCoeffs] = trajectoryGenerator(pWpts, vWpts, aWpts, t, dims, order, corridorExtremes, l, hyperplanes)
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
     
    if isempty(corridorExtremes)
        %% solve equality constraint convex optimization problem
        %% here I tried two solutions
        tic
%         KKT = [Q Aeq';Aeq zeros(size(Q,1)+size(Aeq,1)-size(Aeq,2))];
%         xsol = inv(KKT)*[zeros(size(Q,1),1);beqs(:,1)];
%         ysol = inv(KKT)*[zeros(size(Q,1),1);beqs(:,2)];
    %     zsol = inv(KKT)*[zeros(size(Q,1),1);beqs(:,3)];
        Q1 = blkdiag(Q,Q);
        Aeq1 = blkdiag(Aeq,Aeq);
        KKT = [Q1 Aeq1';Aeq1 zeros(size(Q1,1)+size(Aeq1,1)-size(Aeq1,2))];
        sols = inv(KKT)*[zeros(size(Q1,1),1);[beqs(:,1);beqs(:,2)]];
        xsol = sols(1:size(Q,1));
        ysol = sols(size(Q,1)+1:size(Q,1)*2);
%         zsol = sols(size(Q,1)*2+1:size(Q,1)*3);
        toc
    %     tic
    %     xsol = quadprog(Q,[],[],[],Aeq,beqs(:,1));
    %     ysol = quadprog(Q,[],[],[],Aeq,beqs(:,2));
    %     toc
    else
        %% prepare inequality constraints
        Aineq = [];
        bineq = [];
        for i = 1:size(corridorExtremes,1)
            hyperplane = hyperplanes(i,:);
            numIneq = size(corridorExtremes{i},2);
            Ax = zeros(1,numCoeff*segments);
            Ay = zeros(1,numCoeff*segments);
            Az = zeros(1,numCoeff*segments);
            
            for j = 1:numIneq
                ts = corridorExtremes{i}(j);
                Ax(1,(i-1)*numCoeff+1:i*numCoeff) = hyperplane(1).*[1 ts ts^2 ts^3 ts^4 ts^5];
                Ay(1,(i-1)*numCoeff+1:i*numCoeff) = hyperplane(2).*[1 ts ts^2 ts^3 ts^4 ts^5];
                Az(1,(i-1)*numCoeff+1:i*numCoeff) = hyperplane(3).*[1 ts ts^2 ts^3 ts^4 ts^5];
                %% 1st
                Aineq = [Aineq;[Ax Ay Az];[-Ax -Ay -Az]];
                bineq = [bineq;l-hyperplane(4);l+hyperplane(4)];
                
                Ax(1,(i-1)*numCoeff+1:i*numCoeff) = hyperplane(5).*[1 ts ts^2 ts^3 ts^4 ts^5];
                Ay(1,(i-1)*numCoeff+1:i*numCoeff) = hyperplane(6).*[1 ts ts^2 ts^3 ts^4 ts^5];
                Az(1,(i-1)*numCoeff+1:i*numCoeff) = hyperplane(7).*[1 ts ts^2 ts^3 ts^4 ts^5];
                Aineq = [Aineq;[Ax Ay Az];[-Ax -Ay -Az]];
                bineq = [bineq;l-hyperplane(8);l+hyperplane(8)];
            end
            
            
        end
        Q1 = blkdiag(Q,Q,Q);
        Aeq1 = blkdiag(Aeq,Aeq,Aeq);
        beq1 = [beqs(:,1);beqs(:,2);beqs(:,3)];
        sols = quadprog(2*Q1,[],Aineq,bineq,Aeq1,beq1);
        xsol = sols(1:size(Q,1));
        ysol = sols(size(Q,1)+1:size(Q,1)*2);
        zsol = sols(size(Q,1)*2+1:size(Q,1)*3);
    end

    
%     tic
%     xsol = solveViaQR(Q,Aeq,beqs(:,1));
% %     toc
%     ysol = solveViaQR(Q,Aeq,beqs(:,2));
    polyCoeffs(:,1) = xsol(1:size(Q,1));
    polyCoeffs(:,2) = ysol(1:size(Q,1));
%     polyCoeffs(:,3) = zsol(1:size(Q,1));
end