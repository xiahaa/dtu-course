function ex5_modeling_error_sources
    close all;
    clear all;
    
    % add utils
    addpath('../utils/');
    % sp3 file
    [filename,pathname] = uigetfile('*.sp3','Open a SP3 file...');
    path=strcat(pathname,filename);
    disp(path);
    % parsing sp3 file
    content = sp3fileParser(path);
    % display 
    for i = 1:size(content.sections,1)
        list=sprintf('%d-%d-%d:%d:%d:%f',content.sections{i}.year, ...
            content.sections{i}.Month, content.sections{i}.Day, ...
            content.sections{i}.Hour, content.sections{i}.Minute, ...
            content.sections{i}.Second);
        disp(list);
    end

    % given epoch
    epoch.year = 2018;
    epoch.Month = 9;
    epoch.Day = 17;
    epoch.Hour = 12;
    epoch.Minute = 0;
    epoch.Second = 5;
    
    %% read tpday's information
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
    
    %% parse content
    satInfo = cell(content.satNum,1);
    for i = 1:size(dataOfToday,1)
        data = dataOfToday{i,1};
        time = data.Hour * 3600 + data.Minute * 60 + data.Second;
        for j = 1:size(data.satPos,1)
            satinfo = [time, data.satPos(j).x, data.satPos(j).y, data.satPos(j).z,data.satPos(j).clock];
            satInfo{j}=[satInfo{j};satinfo];
        end
    end 
    
    %%
    ts = epoch.Hour * 3600 + epoch.Minute * 60 + epoch.Second;
    id = 2;
    d_satclk = sat_clock_error(satInfo, id, ts);
    
end