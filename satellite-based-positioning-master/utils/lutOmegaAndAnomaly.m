function [OmegaDeg, AnomalyDeg] = lutOmegaAndAnomaly(slotID)
%lutOmegaAndAnomaly: lookup ascension and anomaly
%   [OmegaDeg, AnomalyDeg] = lutOmegaAndAnomaly(slotID) returns 
%   ascension and anomaly for a given satellite.
%   slotID: satellite ID, now support all 24 satellites of GPS.
%   OmegaDeg: ascension in degree.
%   AnomalyDeg: anomaly in degree.
%   Author: xiahaa@space.dtu.dk
    i = 0; j = 0;
    if slotID(1) == 'A' || slotID(1) == 'a'
        i = 1;
    elseif slotID(1) == 'B' || slotID(1) == 'b'
        i = 2;
    elseif slotID(1) == 'C' || slotID(1) == 'c'
        i = 3;
    elseif slotID(1) == 'D' || slotID(1) == 'd'
        i = 4;
    elseif slotID(1) == 'E' || slotID(1) == 'e'
        i = 5;
    elseif slotID(1) == 'F' || slotID(1) == 'f'
        i = 6;
    end
    j = str2num(slotID(2));
    
    %% lookup table
    lutOmega = [272.85, 332.85, 32.85, 92.85, 152.85, 212.85];
    lutAnomaly = [11.68 41.81 161.79 268.13; ...
                  80.96 173.34 204.38 309.98; ...
                  111.88 241.57 339.67 11.80; ...
                  135.27 167.36 265.45 35.16; ...
                  197.05 302.60 333.69 66.07; ...
                  238.89 345.23 105.21 135.35];
    
    OmegaDeg = lutOmega(i);
    AnomalyDeg = lutAnomaly(i,j);
end