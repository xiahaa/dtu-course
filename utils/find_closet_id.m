function closet_id = find_closet_id(dataOfToday, epoch)
    clockErr = [];
    for i = 1:size(dataOfToday,1)
        clockErr(i) = (dataOfToday{i}.Hour - epoch.Hour)*3600 + ...
                      (dataOfToday{i}.Minute - epoch.Minute)*60 + ...
                      (dataOfToday{i}.Second - epoch.Second);
    end    
    [minval, minid] = min(abs(clockErr));
    closet_id = minid;
end