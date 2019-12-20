function dataOfToday = find_data(content, epoch)
% find today's data
    dataOfToday = {};
    j = 1;
    for i = 1:size(content.sections,1)
        timediff = abs(content.sections{i}.year - epoch.year) * 365 + ...
                   abs(content.sections{i}.Month - epoch.Month) * 30 + ...
                   abs(content.sections{i}.Day - epoch.Day);
        % not today, continue
        if timediff >= 1
            continue;
        end
        % today's data
        da.Hour = content.sections{i}.Hour;
        da.Minute = content.sections{i}.Minute;
        da.Second = content.sections{i}.Second;
        da.satPos = content.sections{i}.satPos;
        dataOfToday{j,1} = da;
        j = j + 1;
    end
end