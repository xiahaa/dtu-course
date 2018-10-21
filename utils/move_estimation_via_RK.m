function s = move_estimation_via_RK(slotID, t0, t)
%move_estimation_via_RK: estimate satellite movement using Rouge-Kutta
%integration.
%   slotID: satellite ID, now support all 24 satellites of GPS.
%   t0: initial time, 0.
%   t: time interval starting from t0.
%   s: return value, estimated movements in meters.
%   Author: xiahaa@space.dtu.dk
    format long;
    
    a = 26559.8*1e3;% semi-major axis
    e = 0;% eccentricity
    inc = deg2rad(55); % GPS inclination
    argPerigee = deg2rad(0);% argument perigee
    % look-up corresponding information for given satellite
    [OmegaDeg, AnomalyDeg] = lutOmegaAndAnomaly(slotID);
    % conversion
    Omega = deg2rad(OmegaDeg);
    M0 = deg2rad(AnomalyDeg);
    % earth gravitational constant
    GM = 3986004.418*1e8;
    
    Kep   = [a, e, inc, Omega, argPerigee, M0]';    % (a,e,i,Omega,omega,M)
    dt = t0;
    y0 = state ( GM, Kep, dt );
    
    % step size
    h = 0.002; % [s]
    % initial values
    te = (t-t0); % end time [s]
    steps = te/h;
    s = 0;

    for i = 1:steps
        y = RK4(@deriv, t0, y0, h, GM);
        y0 = y;
        t0 = t0+h;
        s = s+norm(y(4:6))*h;   
    end
    
end
  

function y = state ( GM, Kep, dt )
    % Keplerian elements at epoch  
    a = Kep(1);  Omega = Kep(4);
    e = Kep(2);  omega = Kep(5);
    i = Kep(3);  M0    = Kep(6);
    % Mean anomaly  
    if (dt==0.0)
        M = M0;
    else
        meanMotion = sqrt(GM)*a^(-3/2);% mean motion
        % mean anomaly
        M = M0 + meanMotion*dt;
    end
    
    % accentric anomaly
    E = estimateEccAnomaly(M, e);
    
    cosE = cos(E);
    sinE = sin(E);
    % Perifocal coordinates
    fac = sqrt ( (1.0-e)*(1.0+e) );
    R = a*(1.0-e*cosE);  % Distance
    V = sqrt(GM*a)/R;    % Velocity
    r = [ a*(cosE-e), a*fac*sinE , 0.0 ]';
    v = [ -V*sinE   , +V*fac*cosE, 0.0 ]';
    % Transformation to reference system (Gaussian vectors)
    Rsq = consR(-Omega,3) * consR(-i,1) * consR(-omega,3);
    r = Rsq*r;
    v = Rsq*v;
    y = [r;v];
end

function yp = deriv (t, y, GM)
    % State vector derivative
    r = y(1:3);
    v = y(4:6);
    yp = [v;-r.*GM/((norm(r))^3)];
end

function y = RK4(func, t, y, h, GM)
    k1 = func(t,y,GM);
    k2 = func(t+h/2,y+(h/2)*k1,GM);
    k3 = func(t+h/2,y+(h/2)*k2,GM);
    k4 = func(t+h,y+h*k3,GM);
    y = y+1/6*k1+1/3*k2+1/3*k3+1/6*k4;
end