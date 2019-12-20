function varargout = solveBySteepestDescentMethod(prs, sat_pos, x0, options)
    tol = 1e6;
    maxIter = 1e6;
    iter = 1;
    oldtol = 1e-6;
    x = x0;
    while abs(tol - oldtol) > 1e-10 && iter < maxIter
        f = evalfunc(x,sat_pos,prs);
        jaco = calcJacobian(x, sat_pos);
        grad = jaco'*f;
        h = -grad;
        a = lineSearch(x, h, sat_pos, prs);
        x = x + h * a;
        iter = iter + 1;
        oldtol = tol;
        tol = f'*f;
    end
    
    f = evalfunc(x,sat_pos,prs);
    RSS = f'*f;
    free = size(prs,1)-size(x0,1);
    MSE = RSS / free;%% s02
    RMSE = sqrt(MSE); %% s0
    
    J = calcJacobian(x, sat_pos);
    Ninv = inv(J'*J);
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
    
end

function a = lineSearch(x0, h, pa, r)
    a = 1;
    beta = 0.9;
    f0 = evalfunc(x0,pa,r);
    F0 = f0'*f0;
    k = 1;
    while 1
        xnew = x0 + a.*h;
        f = evalfunc(xnew,pa,r);
        F = f'*f;
        if F < F0 * 0.95
            break;
        else
            a = a * beta;
            k = k + 1;
        end
        if k > 1e5
            a = 0;
            break;
        end
    end
end

function f = evalfunc(x0, xi, prs)
    f = sqrt((xi(:,1)-x0(1)).^2+(xi(:,2)-x0(2)).^2+(xi(:,3)-x0(3)).^2)+x0(4);
    f = prs - f;
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