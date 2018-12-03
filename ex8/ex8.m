function ex8
    close all;
    clc;
    global cspd;
    cspd = 299792458;% m / s;

    %% add utils path
    addpath('../utils/');
    addpath('../utils/3rdparty');
    
    %% parsing RINEX file
    debug = 1;
    if debug == 0
        [filename,pathname] = uigetfile('*.05O','Open a RINEX file...');
        path=strcat(pathname,filename);
        disp(path);
        contentRINEX = RINEX_file_parser(path);
        save('debug.mat','contentRINEX');
    else
        load('debug.mat');
    end
    date.Year = 2005; date.Month = 9; date.Day = 20; ...
    date.Hour = 3; date.Minute = 28; date.Second = 0;
    %% extract pseudoranges
    data = extract_data(contentRINEX, date);
    [prs, satIDs] = extract_pr(data);

    %% parsing se3 file
    if debug == 0
        [filename,pathname] = uigetfile('*.sp3','Open a SP3 file...');
        path=strcat(pathname,filename);
        disp(path);
        % parsing sp3 file
        contentSP3 = sp3fileParser(path);
        % display 
        for i = 1:size(contentSP3.sections,1)
            list=sprintf('%d-%d-%d:%d:%d:%f',contentSP3.sections{i}.year, ...
                contentSP3.sections{i}.Month, contentSP3.sections{i}.Day, ...
                contentSP3.sections{i}.Hour, contentSP3.sections{i}.Minute, ...
                contentSP3.sections{i}.Second);
            disp(list);
        end
        save('debugsp3.mat','contentSP3');
    else
        load('debugsp3.mat');
    end

    %% epoch
    epoch.year = 2005;
    epoch.Month = 9;
    epoch.Day = 20;
    epoch.Hour = 3;
    epoch.Minute = 28;
    epoch.Second = 0;

    dataOfToday = find_data(contentSP3, epoch);
    satInfo = parse_sat_info(dataOfToday, contentSP3.satNum);
    sp3SatIDs = extract_sat_IDs(contentSP3.satNames, contentSP3.satNum);

    %% estimate transmission time and interpolate the sat positions
    cnt = 2;
    recv_time = date.Hour*3600 + date.Minute*60 + date.Second;
    satPosPred = zeros(numel(prs),3);
    satClkPred = zeros(numel(prs),1);
    for i = 1:numel(prs)
        pr = prs(i);
        transTime = recv_time - pr / cspd;% find transmission time
        neighborIDs = find_neighbor_ids_with_t(dataOfToday, recv_time, cnt);% find nearest id
        satId = satIDs(i);% raw id
        corresId = (find(sp3SatIDs == satId,1));% corresponding id
        [sx,sy,sz] = interp_sat_pos(corresId, neighborIDs, dataOfToday, transTime, '1');% interpolate the satellite positions
        satPosPred(i,:) = [sx,sy,sz];
        
        satClkPred(i) = interp_sat_clk_err(satInfo, corresId, transTime);
    end
    
    %% correct for the coodinate rotation
    we = 7.2921151467e-5;% earth rotational angular rate: rad/s
    for i = 1:numel(prs)
        pr = prs(i);
        transTime = pr / cspd;% find transmission time
        radwe = we * transTime;
        Re = consR(radwe, 3);
        satPosPred(i,:) = (Re * satPosPred(i,:)')';
    end
    %% preliminary position
    x0 = contentRINEX.recPosRaw;
    consParams = struct('a',6378137.0,'f',1/298.257223563); % some constants
    [lat,lon,h] = Cartesian2llh(x0(1),x0(2),x0(3),consParams);
    TECU = 10;%TECU
    %% estimate rough zenith, elevation and azimuth
    azimuths = zeros(numel(prs),1);
    zeniths = zeros(numel(prs),1);
    elevations = zeros(numel(prs),1);   
    dIonos = zeros(numel(prs),1);   
    dTrops = zeros(numel(prs),1); 
    dSatClks = zeros(numel(prs),1); 
    for i = 1:numel(prs)
        dx = satPosPred(i,1) - x0(1);
        dy = satPosPred(i,2) - x0(2);
        dz = satPosPred(i,3) - x0(3);
        %% conversion
        [e,n,u] = WGS842ENU(lat, lon, dx, dy, dz);
        %% compute azimuth and zenith
        [azimuth, zenith, elevation] = calcAzimuthZenithElevation(e,n,u);
        azimuths(i) = azimuth;
        zeniths(i) = zenith;
        elevations(i) = elevation;
        %% correct for atomosphric error
        dIonos(i) = iono_delay_first_order_group(elevation, TECU);
        dTrops(i) = tropo_delay_via_saastamoinen_model(lat, h, elevation, '1');
        %% satellite clock error
        clkerr = satClkPred(i);
        dSatClks(i) = cspd * clkerr * 1e-6;
    end
    
    %% correct
    prs_corr = correct_error_excep_recever(prs, dIonos, dTrops, dSatClks);

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
    options.maxiter = 50;% maximum iteration number
    % two thresholds for iteration: smaller, more iteration.
    options.threshold1 = 1e-12;%  
    options.threshold2 = 1e-12;%
    options.x0_prior = [x0 + 100*ones(3,1);0];

    %% solve with fully corrected pseudorange
    sprior2 = 10^2; %5^2; %prior variance [m^2]
    options.prs_var = sprior2;% prior covariance
    [x_cor,std_x_cor,QDOP_cor,Qenu_cor,llh3] = navSolver(prs_corr, satPosPred, options);

    err_cor = x_cor(1:3) - x0;
    PDOP_cor = sqrt(trace(QDOP_cor(1:3,1:3)));
    TDOP_cor = sqrt(QDOP_cor(4,4));
    HDOP_cor = sqrt(trace(Qenu_cor(1:2,1:2)));
    VDOP_cor = sqrt(Qenu_cor(3,3));

    disp('--------------------Fully Corrected---------------------------');
    disp(strcat("RINEX: ", num2str(x0')));
    disp(strcat("estimation:", num2str(x_cor')));
    disp(strcat("error xyz with raw pr:", num2str(err_cor)));
    disp(strcat("error norm with raw pr:", num2str(norm(err_cor))));    
    disp(strcat('std_x_raw: ',num2str(std_x_cor)));
    disp(strcat('PDOP_raw: ',num2str(PDOP_cor)));
    disp(strcat('TDOP_raw: ',num2str(TDOP_cor)));
    disp(strcat('HDOP_raw: ',num2str(HDOP_cor)));
    disp(strcat('VDOP_raw: ',num2str(VDOP_cor)));

    figure
    hsky = skyPlot(azimuths,zeniths.*180./pi,satIDs,'o');
    set(hsky,'LineWidth',2);
    set(hsky,'MarkerEdgeColor','r');
    set(hsky,'MarkerFaceColor','r');
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

function satInfo = parse_sat_info(dataOfToday, satNum)
    satInfo = cell(satNum,1);
    for i = 1:size(dataOfToday,1)
        data = dataOfToday{i,1};
        time = data.Hour * 3600 + data.Minute * 60 + data.Second;
        for j = 1:size(data.satPos,1)
            satinfo = [time, data.satPos(j).x, data.satPos(j).y, data.satPos(j).z,data.satPos(j).clock];
            satInfo{j}=[satInfo{j};satinfo];
        end
    end
end

function satClkErr = interp_sat_clk_err(sp3SatInfo, id, ts)
%% find the corresponding sp3 indices
    tsp3 = sp3SatInfo{id}(:,1);
    clkrec = sp3SatInfo{id}(:,5);
    %% interpolation the clkrec
    satClkErr = interp1(tsp3,clkrec,ts,'linear');
    satClkErr = satClkErr;
end

function sp3SatIDs = extract_sat_IDs(sat_names, sat_num)
    for i = 1:sat_num
        sp3SatIDs(i) = str2num(sat_names(i,2:end));
    end
end

function varargout = extract_pr(data)
%% here, extract the pseudoranges (L1-L2) from RINEX file.
%  todo: whether or not we need to correct for the Doppler shift.
    prs = zeros(size(data.SNRs,1),1);
    satIDs = zeros(size(data.SNRs,1),1);
    for i = 1:size(data.SNRs,1)
        prs(i) = data.satPos(i).C1;
        satIDs(i) = str2num(data.SNRs{i}(2:end));
    end
    varargout{1} = prs;
    varargout{2} = satIDs;
end

function varargout = extract_data(content, date)
%% extract data corresponding to the given date from the raw and complete
%  RINEX content.
% Author: xiahaa@space.dtu.dk
    dy = content.date.Year - date.Year;
    dm = content.date.Month - date.Month;
    dd = content.date.Day - date.Day;
    if abs(dy)>1 || abs(dm)>1 || abs(dd)>1
        error('The date does not correspond to the given RINEX file!!!!');
        varargout{1} = [];
        return;
    end
    ts1 = date.Hour * 3600 + date.Minute * 60 + date.Second;
    nearest = -1;
    nesrestDt = 1e6;
    for i = 1:size(content.sections,1)
        section = content.sections{i};
        ts2 = section.Hour * 3600 + section.Minute * 60 + section.Second;

        if abs(ts2 - ts1) < nesrestDt
            nesrestDt = abs(ts2 - ts1);
            nearest = i;
        end
    end
    
    if nesrestDt < 1
        data = content.sections{nearest};
        varargout{1} = data;
        disp('found!');
    else
        data = content.sections{nearest};
        varargout{1} = data;
        warning('Only found the nearest!');
    end
end
