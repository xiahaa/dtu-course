function ex7_position_estimation
    close all;
    clc;
    %% add utils path
    addpath('../utils/');
    addpath('../utils/3rdparty');
    
    %% generate test case
    [prs, sat_pos, xreal] = generate_test_case1();
    
    % raw pseudorange
    pseudorange_raw = correct_error_excep_recever(prs(:,1), [], [], []);
    % pseudoranges with satellite clock errors corrected
    pseudorange_partial_cor = correct_error_excep_recever(prs(:,1), [], [], prs(:,5));
    % pseudoranges with atomospheric delays and satellite clock error
    % corrected
    pseudorange_full_cor = correct_error_excep_recever(prs(:,1), prs(:,3), prs(:,4), prs(:,5));
    
    % initial value
%     x0 = xreal + 1000;% this is the case required by the assignemnt
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
    
    options.x0_prior = x0;
    %% solve with raw pseudoranges
    sprior2 = 100000^2; %5^2; %prior variance [m^2]
    options.prs_var = sprior2;% prior covariance
    [x_raw1,std_x_raw1,QDOP_raw1,Qenu_raw1,llh1] = navSolver(pseudorange_raw, sat_pos, options);
    %% solve with partial corrected pseudoranges
    sprior2 = 10^2; %5^2; %prior variance [m^2]
    options.prs_var = sprior2;% prior covariance
    [x_raw2,std_x_raw2,QDOP_raw2,Qenu_raw2,llh2] = navSolverAug(pseudorange_partial_cor, sat_pos, options);
    %% solve with fully corrected pseudorange
    sprior2 = 10^2; %5^2; %prior variance [m^2]
    options.prs_var = sprior2;% prior covariance
    [x_cor,std_x_cor,QDOP_cor,Qenu_cor,llh3] = navSolver(pseudorange_full_cor, sat_pos, options);
    
    err_raw1 = x_raw1(1:3) - xreal;
    err_recerr_raw1 = x_raw1(4) - 0.1*1e-3;
    PDOP_raw1 = sqrt(trace(QDOP_raw1(1:3,1:3)));
    TDOP_raw1 = sqrt(QDOP_raw1(4,4));
    HDOP_raw1 = sqrt(trace(Qenu_raw1(1:2,1:2)));
    VDOP_raw1 = sqrt(Qenu_raw1(3,3));
    
    err_raw2 = x_raw2(1:3) - xreal;
    err_recerr_raw2 = x_raw2(4) - 0.1*1e-3;
    PDOP_raw2 = sqrt(trace(QDOP_raw2(1:3,1:3)));
    TDOP_raw2 = sqrt(QDOP_raw2(4,4));
    HDOP_raw2 = sqrt(trace(Qenu_raw2(1:2,1:2)));
    VDOP_raw2 = sqrt(Qenu_raw2(3,3));
    
    err_cor = x_cor(1:3) - xreal;
    err_recerr_cor = x_cor(4) - 0.1*1e-3;
    PDOP_cor = sqrt(trace(QDOP_cor(1:3,1:3)));
    TDOP_cor = sqrt(QDOP_cor(4,4));
    HDOP_cor = sqrt(trace(Qenu_cor(1:2,1:2)));
    VDOP_cor = sqrt(Qenu_cor(3,3));
    
    disp('--------------------RAW---------------------------');
    disp(strcat("truth: ", num2str(xreal')));
    disp(strcat("estimation:", num2str(x_raw1')));
    disp(strcat("error xyz with raw pr:", num2str(err_raw1)));
    disp(strcat("error norm with raw pr:", num2str(norm(err_raw1))));    
    disp(strcat("error receiver clock with raw pr:", num2str(err_recerr_raw1)));
    disp(strcat('std_x_raw: ',num2str(std_x_raw1)));
    disp(strcat('PDOP_raw: ',num2str(PDOP_raw1)));
    disp(strcat('TDOP_raw: ',num2str(TDOP_raw1)));
    disp(strcat('HDOP_raw: ',num2str(HDOP_raw1)));
    disp(strcat('VDOP_raw: ',num2str(VDOP_raw1)));
    
    disp('--------------------Partial Correction But Augmented State---------------------------');
    disp(strcat("truth: ", num2str(xreal')));
    disp(strcat("estimation:", num2str(x_raw2')));
    disp(strcat("error xyz with raw pr:", num2str(err_raw2)));
    disp(strcat("error norm with raw pr:", num2str(norm(err_raw2))));    
    disp(strcat("error receiver clock with raw pr:", num2str(err_recerr_raw2)));
    disp(strcat('std_x_raw: ',num2str(std_x_raw2)));
    disp(strcat('PDOP_raw: ',num2str(PDOP_raw2)));
    disp(strcat('TDOP_raw: ',num2str(TDOP_raw2)));
    disp(strcat('HDOP_raw: ',num2str(HDOP_raw2)));
    disp(strcat('VDOP_raw: ',num2str(VDOP_raw2)));

    disp('--------------------Fully Corrected---------------------------');
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
    
    disp('--------------------Ionospheric&Tropospheric Delay Estimated---------------------------');
    dx = sat_pos(:,1) - x_cor(1);
    dy = sat_pos(:,2) - x_cor(2);
    dz = sat_pos(:,3) - x_cor(3);
    consParams = struct('a',6378137.0,'f',1/298.257223563); % some constants
    s1 = [];
    s2 = [];
    RE = 6371e3;
    hI = 350e3;
    [lat, lon, height] = Cartesian2llh(x_cor(1),x_cor(2),x_cor(3),consParams);   
    azs = [];
    zns = [];
    PRN = [];
    for i = 1:size(sat_pos,1)
        %% conversion
        [e,n,u] = WGS842ENU(lat, lon, dx(i), dy(i), dz(i));
        %% compute azimuth and zenith
        [azimuth, zenith, elevation] = calcAzimuthZenithElevation(e,n,u);
        zenith = pi/2 - deg2rad(elevation);
        OF = (1-((RE*sin(zenith))/(RE+hI))^2)^(-1/2);
        s1 = [s1;OF];
        s2 = [s2;1/sin(deg2rad(elevation))];
        
        azs = [azs;azimuth];
        zns = [zns;zenith];
        PRN = [PRN;i];
    end
    A = [s1 s2];
    baug = A * [x_raw2(5);x_raw2(6)];
    disp('Truth');
    disp([prs(:,3), prs(:,4)]);
    disp('Estimated');
    disp(baug);
    
    figure
    hsky = skyPlot(azs,zns.*180./pi,PRN,'o');
    set(hsky,'LineWidth',2);
    set(hsky,'MarkerEdgeColor','r');
    set(hsky,'MarkerFaceColor','r')
%     figure
    % sky plot
%     polar(azs(:),zns(:).*180./pi,'ob');
    %h = mmpolar(azim(I),zen_ang(I).*180/pi,'*b','TZeroDirection','North','RLimit',[0 90]);
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

function [prs, sat_pos, x0] = generate_test_case1()
%% generate_test_case1: 
%   generate pseudoranges for a hard-coded postision (DTU 101).
%   inputs: none
%   outputs:
%       prs: pseudoranges
%       sat_pos: satellite positions
%       x0: true receiver position
%% Author: xiahaa@space.dtu.dk

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
