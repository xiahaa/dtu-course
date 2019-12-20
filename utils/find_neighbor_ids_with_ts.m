function neighbor_ids = find_neighbor_ids_with_ts(dataOfToday, ts, cnt)
    clockErr = [];
    for i = 1:size(dataOfToday,1)
        clockErr(i) = (dataOfToday{i}.Hour*3600 + dataOfToday{i}.Minute*60 + dataOfToday{i}.Second) - ts;
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