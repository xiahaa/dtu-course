function utc = utcdate(jd)
%% convert from julian date to utc
% Author: xiahaa@space.dtu.dk
    z = fix(jd + .5);
    fday = jd + .5 - z;
    if (fday < 0)
        fday = fday + 1;
        z = z - 1;
    end
    if (z < 2299161)
        a = z;
    else
        alpha = floor((z - 1867216.25) / 36524.25);
        a = z + 1 + alpha - floor(alpha / 4);
    end

    b = a + 1524;
    c = fix((b - 122.1) / 365.25);
    d = fix(365.25 * c);
    e = fix((b - d) / 30.6001);
    utc.day = b - d - fix(30.6001 * e) + fday;

    if (e < 14)
       utc.month = e - 1;
    else
       utc.month = e - 13;
    end
 
    if (utc.month > 2)
       utc.year = c - 4716;
    else
       utc.year = c - 4715;
    end
    utc.hour = abs(utc.day-floor(utc.day))*24;
    utc.minute = abs(utc.hour-floor(utc.hour))*60;
    utc.second = round(abs(utc.minute-floor(utc.minute))*60);
    utc.day = floor(utc.day);
    utc.hour = floor(utc.hour);
    utc.minute = floor(utc.minute);
end