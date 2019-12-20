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
         
    mpol t p0 p1 p2 p3 p4 p5 t1 p10 p11 p12 p13 p14 p15
    
    cost = t^3*(192*p4^2 + 240*p3*p5) + 720*p5^2*t^5 + 36*p3^2*t + 144*p3*p4*t^2 + 720*p4*p5*t^4 + ...
           t1^3*(192*p14^2 + 240*p13*p15) + 720*p15^2*t1^5 + 36*p13^2*t1 + 144*p13*p14*t1^2 + 720*p14*p15*t1^4;

    order = 5;
    %% order of trajectory 5
    numCoeff = order + 1;
    segments = size(pWpts,1)-1;

    %% equality constraint
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
    %% continuity constraints
    for i=2:size(pWpts,1)-1
        if vWpts(i,1) == inf
            alreadySetCons = alreadySetCons + 1;
        end
        if aWpts(i,1) == inf
            alreadySetCons = alreadySetCons + 1;
        end
    end
    
    K = [];
    %% P1
    K1 = p0 == pWpts(1,1);
    K2 = p1 == vWpts(1,1);
    K3 = 2 * p2 == aWpts(1,1);
%     
    K4 = p0 + t * p1 + t^2 * p2 + t^3 * p3 + t^4 * p4 + t^5 * p5 == pWpts(2,1);
%     K5 = p10 == pWpts(2,1);
%     
%     K5 = p10 + t1 * p11 + t1^2 * p12 + t1^3 * p13 + t1^4 * p14 + t1^5 * p15 == pWpts(3,1);
%     K6 =              p11 + 2*t1 * p12 + 3*t1^2 * p13 + 4*t1^3 * p14 + 5*t1^4 * p15 == vWpts(3,1);
%     K7 =                      2 * p12 + 6*t1 * p13 + 12*t1^2 * p14 + 20*t1^3 * p15 == aWpts(3,1);
    
    
%     K8 = p1 + 2*t * p2 + 3*t^2 * p3 + 4*t^3 * p4 + 5*t^4 * p5 - p11 == 0;
%     K9 = 2 * p2 + 6*t * p3 + 12*t^2 * p4 + 20*t^3 * p5 - 2 * p12 == 0;
    
    K = [K1 K2 K3 K4 t>=1];%K5 K6 K7
    
    Prob = msdp(min(cost),K);
    [status,obj] = msol(Prob);
    xsol = meas;
    
    ts = [0 double(t0) double(t1)];
    polyCoeffs = [double(p0)';double(p1)'];
    
    
    topt = ts;
    
%     [polyCoeffsNew] = trajectoryGenerator(pWpts, vWpts, aWpts, topt, dims, order);
    visulization(pWpts, polyCoeffs, topt, order, vmax, amax);
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


