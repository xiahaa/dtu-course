function [q1, q2, q3, i1, i2, i3] = calcSatPosition(slotID, t0, t)
%calcSatPosition: compute satellite position in orbital and CIRS coordinate frame
%   [q1, q2, q3, i1, i2, i3] = calcSatPosition(slotID, t0, t) compute 
%   satellite position in orbital and CIRS coordinate frame.
%   slotID: satellite ID, now support all 24 satellites of GPS.
%   t0: initial time, 0.
%   t: time interval starting from t0.
%   q1, q2, q3: return value, satellite positions in orbital frame.
%   i1, i2, i3: return value, satellite positions in ECEF frame.
%   Author: xiahaa@space.dtu.dk
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
    meanMotion = sqrt(GM)*a^(-3/2);% mean motion
    % mean anomaly
    M = M0 + meanMotion*(t-t0);
    % accentric anomaly
    E = estimateEccAnomaly(M, e);
    % orbital positions
    q1 = a*cos(E) - a*e;
    q2 = a*sqrt(1-e^2)*sin(E);
    q3 = 0;
    
    % constuct rotation matrix
    ROmega = consR(Omega,3);
    Rinc = consR(inc,1);
    RargPerigee = consR(argPerigee,3);
    
%     ROmega = [cos(Omega) sin(Omega) 0;
%               -sin(Omega) cos(Omega) 0;
%               0 0 1];
%     Rinc = [1 0 0;
%             0 cos(inc) sin(inc);
%             0 -sin(inc) cos(inc)];
%     RargPerigee = [cos(argPerigee) sin(argPerigee) 0;
%               -sin(argPerigee) cos(argPerigee) 0;
%               0 0 1];

    Rqs = RargPerigee * Rinc * ROmega;
    Rsq = Rqs';
    % transformation to ECEF
    v = Rsq * [q1;q2;q3];
    i1 = v(1);
    i2 = v(2);
    i3 = v(3);
end