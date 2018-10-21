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

    % local position llh
    testCase(1,:) = [55.78575300466123,12.525384183973078,40];% DTU 101
    
    % given epoch
    epoch.year = 2018;
    epoch.Month = 9;
    epoch.Day = 17;
    epoch.Hour = 12;
    epoch.Minute = 0;
    epoch.Second = 5;
    
    %% read tpday's information
    dataOfToday = find_data(content, epoch);   
    
    %% satellite position prediction
    satPosPred = zeros(content.satNum,3);
    t_interp = epoch.Hour*3600+epoch.Minute*60+ epoch.Second;
    cnt = 2;
    neighbor_ids = find_neighbor_ids(dataOfToday, epoch, cnt);
    for i = 1:content.satNum
        satId = i;
        [sx,sy,sz] = interp_sat_pos(satId, neighbor_ids, dataOfToday, t_interp, '1');
        satPosPred(i,:) = [sx,sy,sz];
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
    ts = t_interp;
    TECU = 10;%TECU
    recerr = 0.1;%ms
    prs = [];
    for i = 1:content.satNum
        [visibility, pr, R, d_iono, d_trop, d_satclk, d_recclk] = calc_pseudorange(testCase(1,1), testCase(1,2), testCase(1,3), ...
                                      satPosPred(i,:), satInfo, i, ts, TECU, recerr);
        if visibility == 1
            prs = [prs;[pr R, d_iono, d_trop, d_satclk, d_recclk]];
        end
    end
    
    disp(prs);
    
end