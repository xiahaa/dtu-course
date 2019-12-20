function E = estimateEccAnomaly(M, e)
%estimateEccAnomaly: estimate eccentric anomaly iteratively.
%   E = estimateEccAnomaly(M, e) returns estimated eccentric anomaly.
%   M: mean anomaly.
%   e: eccentricity.
%   E: eccentric anomaly.
%   Author:xiahaa@space.dtu.dk
    E = 0;
    iter = 1;
    maxIter = 1e3;
    error = 1e6;
    while error > 1e-3 && iter <= maxIter
        newE = M + e*sin(E);
        error = abs(newE) - E;
        E = newE;
    end
end