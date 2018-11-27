function ex7_dop_evaluation
    close all;
    clear all;
    clc;
    %% add utils path
    addpath('../utils/');
    addpath('../utils/3rdparty');
    
    [rec_llh, rec_xyz, content] = prepare();
    x0 = [0;0;0];
    x0 = [x0;0];%% add initial receiver clock error
    
    %% options for the navSolver
    %   for initialization:
    %       1: use the x0 provided by the user;
    %       2: use the proposed second-order-cone programming;
    %       3: use the proposed direct linear transformation method;
    %       4: use the Bancraft method;
    %   for iterative solver:
    %       1: Gauss-Newton;
    %       2: Steepest Dscent;
    %       3: Levenberg-Marquardt;   
    options.initialization = 1;% can be 1,2,3,4
    options.solver = 1;% can be 1,2,3
    options.verbose = 0;% will generate intermediate print info
    options.maxiter = 100;% maximum iteration number
    %% two thresholds for iteration: smaller, more iteration.
    options.threshold1 = 1e-6;%  
    options.threshold2 = 1e-6;%
    
    % initial epoch
    epoch.year = 2018;
    epoch.Month = 9;
    epoch.Day = 16;
    epoch.Hour = 12;
    epoch.Minute = 0;
    epoch.Second = 0;
    
    PDOPs = [];
    HDOPs = [];
    GDOPs = [];
    VDOPs = [];
    kk = 1;
    times = [];
    for t = 0:60*15:3600*24
        cepoch.year = 2018;
        cepoch.Month = 9;
        cepoch.Day = epoch.Day + floor((epoch.Hour + floor(t/(3600)))/24);
        cepoch.Hour = mod(epoch.Hour + floor(t/(3600)),24);
        cepoch.Minute = mod(epoch.Minute + floor(t/(60)),60);
        cepoch.Second = epoch.Second + mod(t,60);
%         disp(cepoch);
        %% generate test case
        [prs, sat_pos] = find_sat_xyz(cepoch,rec_llh, content);
        % pseudoranges with atomospheric delays and satellite clock error
        % corrected
        pseudorange_full_cor = correct_error_excep_recever(prs(:,1), prs(:,3), prs(:,4), prs(:,5));

        options.x0_prior = x0;
        %% solve with fully corrected pseudorange
        sprior2 = 10^2; %5^2; %prior variance [m^2]
        options.prs_var = sprior2;% prior covariance
        [x_cor,std_x_cor,QDOP_cor,Qenu_cor,llh3] = navSolver(pseudorange_full_cor, sat_pos, options);
        PDOP_cor = sqrt(trace(QDOP_cor(1:3,1:3)));
        GDOP_cor = sqrt(trace(QDOP_cor(1:4,1:4)));
        HDOP_cor = sqrt(trace(Qenu_cor(1:2,1:2)));
        VDOP_cor = sqrt(trace(Qenu_cor(3,3)));
        
        PDOPs = [PDOPs;PDOP_cor];
        HDOPs = [HDOPs;HDOP_cor];
        GDOPs = [GDOPs;GDOP_cor];
        VDOPs = [VDOPs;VDOP_cor];
        
        sat_pos_all{kk} = sat_pos;
        kk = kk + 1;
        times = [times;t];
    end
        
    figure
    plot(times,PDOPs,'b-o','LineWidth',3);grid on;
    title('PDOP over 24 hours','Interpreter','latex');
    xlabel('time: (s)','Interpreter','latex');
    ylabel('PDOP','Interpreter','latex');

    figure
    plot(times,GDOPs,'b-o','LineWidth',3);grid on;
    title('GDOP over 24 hours','Interpreter','latex');
    xlabel('time: (s)','Interpreter','latex');
    ylabel('GDOP','Interpreter','latex');
    
    figure
    plot(times,HDOPs,'b-o','LineWidth',3);grid on;
    title('HDOP over 24 hours','Interpreter','latex');
    xlabel('time: (s)','Interpreter','latex');
    ylabel('HDOP','Interpreter','latex');
    
    figure
    plot(times,VDOPs,'b-o','LineWidth',3);grid on;
    title('VDOP over 24 hours','Interpreter','latex');
    xlabel('time: (s)','Interpreter','latex');
    ylabel('VDOP','Interpreter','latex');
    
    [val1,minid] = min(GDOPs);
    [val2,maxid] = max(GDOPs);
    
    hmin = drawSkyPlot(sat_pos_all{minid}, rec_xyz);
    title(strcat('Skyplot with the minimum GDOP: ', num2str(val1)),'Interpreter','latex');
    
    hmax = drawSkyPlot(sat_pos_all{maxid}, rec_xyz);
    title(strcat('Skyplot with the maximum GDOP: ', num2str(val2)),'Interpreter','latex');
    
end

function h1 = drawSkyPlot(sat_pos, xyz)
    dx = sat_pos(:,1) - xyz(1);
    dy = sat_pos(:,2) - xyz(2);
    dz = sat_pos(:,3) - xyz(3);
    consParams = struct('a',6378137.0,'f',1/298.257223563); % some constants
    RE = 6371e3;
    hI = 350e3;
    [lat, lon, height] = Cartesian2llh(xyz(1),xyz(2),xyz(3),consParams);   
    azs = [];
    zns = [];
    PRN = [];
    for i = 1:size(sat_pos,1)
        %% conversion
        [e,n,u] = WGS842ENU(lat, lon, dx(i), dy(i), dz(i));
        %% compute azimuth and zenith
        [azimuth, zenith, elevation] = calcAzimuthZenithElevation(e,n,u);
        zenith = pi/2 - deg2rad(elevation);
        
        azs = [azs;azimuth];
        zns = [zns;zenith];
        PRN = [PRN;i];
    end
    
    h1 = figure
    hsky = skyPlot(azs,zns.*180./pi,PRN,'o');
    set(hsky,'LineWidth',2);
    set(hsky,'MarkerEdgeColor','r');
    set(hsky,'MarkerFaceColor','r')
end

function prs = correct_error_excep_recever(prs_raw, d_iono, d_trop, d_satclk)
%% correct_error_excep_recever: 
%   correct pseudoranges with estimated atmospeheric effects.
%   inputs: 
%       prs_raw: raw pseudoranges;
%       d_iono: ionospheric delay;
%       d_trop: tropospheric delay;
%       d_satclk: satellite clock error;
%   outputs:
%       prs: corrected pseudoranges
%% Author: xiahaa@space.dtu.dk

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

function [rec_llh, rec_xyz, content] = prepare()
    % local position llh
    testCase(1,:) = [55.78575300466123,12.525384183973078,0];% DTU 101
    testid = 1;
    consParams = struct('a',6378137.0,'f',1/298.257223563); % some constants
    lat = testCase(testid,1);%
    lon = testCase(testid,2);%
    height = testCase(testid,3);%
    [xo,yo,zo] = llhtoCartesian(lat, lon, height, consParams);% to ECEF
    
    rec_llh = testCase(1,:);
    rec_xyz = [xo,yo,zo];
    
    % sp3 file
    [filename,pathname] = uigetfile('*.sp3','Open a SP3 file...');
    path=strcat(pathname,filename);
    disp(path);
    % parsing sp3 file
    content = sp3fileParser(path);
end

function [prs, sat_pos] = find_sat_xyz(epoch,rec_llh, content)
%% generate_test_case1: 
%   generate pseudoranges for a hard-coded postision (DTU 101).
%   inputs: none
%   outputs:
%       prs: pseudoranges
%       sat_pos: satellite positions
%       x0: true receiver position
%% Author: xiahaa@space.dtu.dk

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
        [visibility, pr, R, d_iono, d_trop, d_satclk, d_recclk, clkerr] = calc_pseudorange(rec_llh(1,1), rec_llh(1,2), rec_llh(1,3), ...
                                      satPosPred(i,:), satInfo, i, ts, TECU, recerr);
        if visibility == 1
            ids = [ids;i];
            prs = [prs;[pr, R, d_iono, d_trop, d_satclk, d_recclk, clkerr]];
        end
    end
    
%     disp(prs);
    
    %% assign output
    sat_pos = satPosPred(ids,:);
end
