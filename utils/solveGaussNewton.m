function varargout = solveGaussNewton(prs, sat_pos, x0, options)        
    %% weight
    x = x0;
    n = size(prs,1);
    p = size(x0,1);
    free = n - p;
    if options.verbose == 1
        %% save some intermediate results
        iters = [];%% iterations cnt
        X = [];%% solution of each iteration
        vtPv = [];%% residual cost
    end
    
    for iter = 1:1:options.maxiter
        %% f(x0)
        f = evalfunc(sat_pos, x);
        %% residuals
        y = prs - f;
        %% jacobian
        J = calcJacobian(x, sat_pos);
        %% A and b
        A = J'*J;
        b = J'*y;
        dxhat = lssolver(A,b);
        %% update
        x = x - dxhat;
        %% verbose
        if options.verbose == 1
            %% intermediate results for residual analysis
            Ninv = inv(A);
            vhat = y;%%residual;
            vtPv = [vtPv vhat'*vhat];%% cost
            X = [X x];
            iters = [iters iter];
        else
            Ninv = inv(A);
            vhat = y;%% new residual;
        end
        
        %% iteration stop criteria
        cond1 = (vhat'*vhat);
        cond2 = norm(dxhat);
        if cond1 < options.threshold1 || cond2 < options.threshold2
            break;
        end
    end
    %% RSS
    RSS = vhat'*vhat;
    MSE = RSS / free;%% s02
    RMSE = sqrt(MSE); %% s0
    
    Dx = MSE * Ninv;
    std_x = sqrt(diag(Dx));
    
    QDOP = computeDOP(J,[],options.prs_var);
    PDOP = sqrt(trace(QDOP(1:3,1:3)));
    TDOP = sqrt(QDOP(4,4));
    GDOP = sqrt(trace(QDOP));
    
    disp(strcat('iter: ', num2str(iter)));
    disp(strcat('xsol: ', num2str(x)));
    disp(strcat('std dev of x: ', num2str(std_x)));
    disp(strcat('PDOP: ', num2str(PDOP)));
    disp(strcat('TDOP: ', num2str(TDOP)));
    disp(strcat('GDOP: ', num2str(GDOP)));
    
    consParams = struct('a',6378137.0,'f',1/298.257223563); % some constants
    [lat, lon, height] = Cartesian2llh(x(1),x(2),x(3),consParams);
    R1 = consR(deg2rad(90-lat),1);
    R3 = consR(deg2rad(90+lon),3);
    R = R1*R3;
    Qenu = R*QDOP(1:3,1:3)*R';
    HDOP = sqrt(trace(Qenu(1:2,1:2)));
    VDOP = sqrt(Qenu(3,3));
    
    disp(strcat('HDOP: ', num2str(HDOP)));
    disp(strcat('VDOP: ', num2str(VDOP)));
    
    varargout{1} = x;
    varargout{2} = std_x;
    varargout{3} = QDOP;
    varargout{4} = Qenu;
    varargout{5} = [lat, lon, height];
    
    if options.verbose == 1
        %% residual analysis    
        %% variance 
        D0 = MSE;
        std_y = sqrt(diag(D0));
        
        Dyhat = MSE * J*Ninv*J';
        std_yhat = sqrt(diag(Dyhat));
        
        Dehat = D0 - Dyhat;
        std_ehat = sqrt(diag(Dehat));
        
        %% compute the correlation
        aux = diag(1./std_x);
        corr_xhat = aux * Dx * aux;
        
        hvpv = display_customized(iters,vtPv,'$v^TPv$',2);
        for i = 1:p
            h(i) = display_customized(iters,X(i,:),strcat('x(', num2str(i),') v.s. iteration'),2);
        end
        
        t = x./std_x;
        pt = betainc(free./(free+t.^2),0.05*free,0.05);
        
        pchi2 = 1-gammainc(0.5*vtPv(end),0.5*free);
        
        disp('probalibity of find large t');
        disp(pt);
        
        disp('probalibity of find large vtPvt');
        disp(pchi2);
        
        %% ellipsoid drawing
        metrics = {'m','m','m'};
        he1 = drawConfidenceEllipsoid(Dx(1:3,1:3), 'chi', metrics);
        he2 = drawConfidenceEllipsoid(Dx(1:3,1:3), 'F', metrics);
        
        he1 = drawConfidenceEllipsoid(Dx(1:3,1:3), 'chi', metrics);
        he2 = drawConfidenceEllipsoid(Dx(1:3,1:3), 'F', metrics);
    end
end

function f = evalfunc(xi, x0)
    f = sqrt((xi(:,1)-x0(1)).^2+(xi(:,2)-x0(2)).^2+(xi(:,3)-x0(3)).^2)+x0(4);
end

function J = calcJacobian(x0, xi)
    r = sqrt((xi(:,1)-x0(1)).^2+(xi(:,2)-x0(2)).^2+(xi(:,3)-x0(3)).^2);
    
    j1x = @(x, xi)((x0(1)-xi(:,1))./r);
    j1y = @(x, xi)((x0(2)-xi(:,2))./r);
    j1z = @(x, xi)((x0(3)-xi(:,3))./r);
    j1dt = @(x, xi) (ones(size(xi,1),1));
    
    n = size(xi,1);
    J = zeros(n, 4);
    
    j11 = j1x(x0,xi);
    j12 = j1y(x0,xi);
    j13 = j1z(x0,xi);
    j14 = j1dt(x0,xi);

    J = [j11 j12 j13 j14];
    J = -J;
end