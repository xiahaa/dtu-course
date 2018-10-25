function ex7_position_estimation
    close all;
    addpath('../utils/');
    addpath('../utils/3rdparty');
    
    [prs, sat_pos, xreal] = generate_test_case1();
    pseudorange_raw = correct_error_excep_recever(prs(:,1), [], [], prs(:,5));
    pseudorange_cor = correct_error_excep_recever(prs(:,1), prs(:,3), prs(:,4), prs(:,5));
    
    x0 = xreal + -100000*rand(3,1)+200000*2;
    x0 = [x0;0];
    
    sprior2 = 10^2; %5^2; %prior variance [m^2]
    options.usePrior = 1;options.useSOCP = 0;options.useDLT = 0;
    options.useWLS = 0;options.useGN = 1;options.useSD = 0;options.useLM = 0;
    options.verbose = 0;
    options.maxiter = 100;
    options.threshold1 = 1e-6;
    options.threshold2 = 1e-6;
    options.prs_var = sprior2;
    
    options.x0_prior = x0;
%     [x_raw,std_x_raw,QDOP_raw,Qenu_raw,llh_raw] = navSolver(pseudorange_raw, sat_pos, options);
    [x_raw,std_x_raw,QDOP_raw,Qenu_raw,llh_raw] = navSolverAug(pseudorange_raw, sat_pos, options);
    
    [x_cor,std_x_cor,QDOP_cor,Qenu_cor,llh_cor] = navSolver(pseudorange_cor, sat_pos, options);
    
    err_raw = x_raw(1:3) - xreal;
    err_recerr_raw = x_raw(4) - 0.1*1e-3;
    PDOP_raw = sqrt(trace(QDOP_raw(1:3,1:3)));
    TDOP_raw = sqrt(QDOP_raw(4,4));
    HDOP_raw = sqrt(trace(Qenu_raw(1:2,1:2)));
    VDOP_raw = sqrt(Qenu_raw(3,3));
    
    err_cor = x_cor(1:3) - xreal;
    err_recerr_cor = x_cor(4) - 0.1*1e-3;
    PDOP_cor = sqrt(trace(QDOP_cor(1:3,1:3)));
    TDOP_cor = sqrt(QDOP_cor(4,4));
    HDOP_cor = sqrt(trace(Qenu_cor(1:2,1:2)));
    VDOP_cor = sqrt(Qenu_cor(3,3));
    
    disp('--------------------RAW---------------------------');
    disp(strcat("truth: ", num2str(xreal')));
    disp(strcat("estimation:", num2str(x_raw')));

    disp(strcat("error xyz with raw pr:", num2str(err_raw)));
    disp(strcat("error norm with raw pr:", num2str(norm(err_raw))));    
    disp(strcat("error receiver clock with raw pr:", num2str(err_recerr_raw)));
    disp(strcat('std_x_raw: ',num2str(std_x_raw)));
    disp(strcat('PDOP_raw: ',num2str(PDOP_raw)));
    disp(strcat('TDOP_raw: ',num2str(TDOP_raw)));
    disp(strcat('HDOP_raw: ',num2str(HDOP_raw)));
    disp(strcat('VDOP_raw: ',num2str(VDOP_raw)));

    disp('--------------------COR---------------------------');
    disp(strcat("truth: ", num2str(xreal')));
    disp(strcat("estimation:", num2str(x_cor')));
    
    disp(strcat("error xyz with raw pr:", num2str(err_cor)));
    disp(strcat("error norm with raw pr:", num2str(norm(err_cor))));    
    disp(strcat("error receiver clock with raw pr:", num2str(err_recerr_cor)));
    disp(strcat('std_x_raw: ',num2str(std_x_cor)));
    disp(strcat('PDOP_raw: ',num2str(PDOP_cor)));
    disp(strcat('TDOP_raw: ',num2str(TDOP_cor)));
    disp(strcat('HDOP_raw: ',num2str(HDOP_cor)));
    disp(strcat('VDOP_raw: ',num2str(VDOP_cor)));
    
end

function prs = correct_error_excep_recever(prs_raw, d_iono, d_trop, d_satclk)
    prs = prs_raw;
    if ~isempty(d_iono)
        prs = prs - d_iono;
    end
    if ~isempty(d_trop)
        prs = prs - d_trop;
    end
    if ~isempty(d_satclk)
        prs = prs + d_satclk;
    end
end

function [prs, sat_pos, x0] = generate_test_case1()
    % local position llh
    testCase(1,:) = [55.78575300466123,12.525384183973078,0];% DTU 101
    consParams = struct('a',6378137.0,'f',1/298.257223563); % some constants
    testid = 1;
    lat = testCase(testid,1);%
    lon = testCase(testid,2);%
    height = testCase(testid,3);%
    [xo,yo,zo] = llhtoCartesian(lat, lon, height, consParams);% to ECEF
    
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
%         disp(list);
    end
    
    % given epoch
    epoch.year = 2018;
    epoch.Month = 9;
    epoch.Day = 16;
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
    ids = [];
    for i = 1:content.satNum
        [visibility, pr, R, d_iono, d_trop, d_satclk, d_recclk, clkerr] = calc_pseudorange(testCase(1,1), testCase(1,2), testCase(1,3), ...
                                      satPosPred(i,:), satInfo, i, ts, TECU, recerr);
        if visibility == 1
            ids = [ids;i];
            prs = [prs;[pr, R, d_iono, d_trop, d_satclk, d_recclk, clkerr]];
        end
    end
    
%     disp(prs);
    
    %% assign output
    x0 = [xo, yo, zo]';
    sat_pos = satPosPred(ids,:);
end
