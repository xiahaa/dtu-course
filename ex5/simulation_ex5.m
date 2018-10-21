function [satPosPred, visibilities, prs] = simulation_ex5(epoch, content, llh, TECU, recerr)
    % add utils
    addpath('../utils/');
    
    %% read tpday's information
    dataOfToday = find_data(content, epoch);   
    
    %% satellite position prediction
    satPosPred = zeros(content.satNum, 3);
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
    prsd = [];prs = [];
    visibilities = zeros(content.satNum,1);
    for i = 1:content.satNum
        [visibility, pr, R, d_iono, d_trop, d_satclk, d_recclk] = ...
            calc_pseudorange(llh(1,1), llh(1,2), llh(1,3), ...
                             satPosPred(i,:), satInfo, i, ts, TECU, recerr);
        if visibility == 1
            visibilities(i) = 1;
            prsd = [prsd;[i, pr, R, d_iono, d_trop, d_satclk, d_recclk]];
        end
        prs = [prs;[i, pr, R, d_iono, d_trop, d_satclk, d_recclk]];
    end
    
    disp(prsd);
    
end