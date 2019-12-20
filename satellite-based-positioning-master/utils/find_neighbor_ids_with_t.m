function neighbor_ids = find_neighbor_ids_with_t(dataOfToday, time, cnt)
    clockErr = [];
    for i = 1:size(dataOfToday,1)
        t1 = dataOfToday{i}.Hour * 3600 + dataOfToday{i}.Minute * 60 + dataOfToday{i}.Second;
        t2 = time;
        clockErr(i) = t1 - t2;
    end    
    [~, minid] = min(abs(clockErr));
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