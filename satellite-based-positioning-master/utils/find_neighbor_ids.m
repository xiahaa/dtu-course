function neighbor_ids = find_neighbor_ids(dataOfToday, epoch, cnt)
    clockErr = [];
    for i = 1:size(dataOfToday,1)
        clockErr(i) = (dataOfToday{i}.Hour - epoch.Hour)*3600 + ...
                      (dataOfToday{i}.Minute - epoch.Minute)*60 + ...
                      (dataOfToday{i}.Second - epoch.Second);
    end    
    [minval, minid] = min(abs(clockErr));
    neighbor_ids = [];
    if (minid - cnt) >= 1 && (minid + cnt) <= numel(clockErr)
        neighbor_ids = (minid - cnt):1:(minid + cnt);
    else
        if (minid - cnt) >= 1
            neighbor_ids = (minid - cnt):1:minid;
        else
            neighbor_ids = minid:1:(minid+cnt);
        end
    end
end