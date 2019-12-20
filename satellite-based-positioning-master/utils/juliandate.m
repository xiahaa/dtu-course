function jd = juliandate(utc)
%% convert from utc to julian date
% Author: xiahaa@space.dtu.dk
    y = utc.year;
    m = utc.month;
    b = 0; c = 0;
    if (m <= 2)
        y = y - 1;
        m = m + 12;
    end
    if (y < 0)
        c = -.75;
    end

    % check for valid calendar date
    if (utc.year < 1582)
       % null
    elseif (utc.year > 1582)
       a = fix(y / 100);
       b = 2 - a + floor(a / 4);
    elseif (utc.month < 10)
       % null
    elseif (utc.month > 10)
       a = fix(y / 100);
       b = 2 - a + floor(a / 4);
    elseif (utc.day <= 4)
       % null
    elseif (utc.day > 14)
       a = fix(y / 100);
       b = 2 - a + floor(a / 4);
    else
        error('\n\n  this is an invalid calendar date!!\n');
    end
    jd = fix(365.25 * y + c) + fix(30.6001 * (m + 1));
    jd = jd + utc.day + b + 1720994.5;
    jd = jd + (utc.hour+utc.minute/60+utc.second/3600)/24;
end