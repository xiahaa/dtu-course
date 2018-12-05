function ex_9_gdop_galileo_gps
    close all;
    clc;clear all;
    %% add utility functions
    addpath('../utils/');
    addpath('../utils/3rdparty');
    
    %% load Galileo Kepler Params
    GalileoKepParams;
    numGalileoSat = size(GalileoKeplerTb,1);    
    
    %% GPS
    satGPS = {'A1','A2','A3','A4', ...
              'B1','B2','B3','B4', ...
              'C1','C2','C3','C4', ...
              'D1','D2','D3','D4', ...
              'E1','E2','E3','E4', ...
              'F1','F2','F3','F4'}';
    numGPSSat = size(satGPS,1);    
    
    %% simulated midnight utc time
    utcMidNight.year = 2018;
    utcMidNight.month = 12;
    utcMidNight.day = 4;
    utcMidNight.hour = 0; 
    utcMidNight.minute = 0;   
    utcMidNight.second = 0;
    
    jdMidNight = juliandate(utcMidNight);
    gast = jd2gast(jdMidNight);%% degree
    
    %% prepare
    [rec_llh, rec_xyz] = prepare();
    lat = rec_llh(1);lon = rec_llh(2);
    
    PDOPs = [];
    times = [];
    for t = 0:60*15:3600*24
        theta = gast + t / 3600 * 15;
        R = consR(theta,3);
        
        GalileosatPos = zeros(numGalileoSat,3);
        GPSsatPos = zeros(numGPSSat,3);
        %% Galileo
        for j = 1:numGalileoSat
            [~, ~, ~, i1, i2, i3] = calc_sat_pos_with_Kepler(GalileoKeplerTb(j,:), 0, t);
            GalileosatPos(j,1) = i1;
            GalileosatPos(j,2) = i2;
            GalileosatPos(j,3) = i3;
        end
        for j = 1:numGPSSat
            slotID = satGPS{j};
            [~, ~, ~, i1, i2, i3] = calcSatPosition(slotID, 0, t);
            GPSsatPos(j,1) = i1;
            GPSsatPos(j,2) = i2;
            GPSsatPos(j,3) = i3;
        end
        %% transform
        GalileosatPosWGS = (R * GalileosatPos')';
        GPSsatPosWGS = (R * GPSsatPos')';
        %% azi, elv
        dGalileo = GalileosatPosWGS - repmat(rec_xyz,numGalileoSat,1);
        dGPS = GPSsatPosWGS - repmat(rec_xyz,numGPSSat,1);
        visibilityGalileo = checkVisibility(numGalileoSat, dGalileo, lat, lon);
        visibilityGPS = checkVisibility(numGPSSat, dGPS, lat, lon);
        %% normalization
        normdGalileo = sqrt(dGalileo(:,1).^2+dGalileo(:,2).^2+dGalileo(:,3).^2);
        normdGPS = sqrt(dGPS(:,1).^2+dGPS(:,2).^2+dGPS(:,3).^2);
        
        idGalileo = (visibilityGalileo == 1);
        idGPS = (visibilityGPS == 1);
        
        ndGalileo = dGalileo(idGalileo,:) ./ normdGalileo(idGalileo);
        ndGPS = dGPS(idGPS,:) ./ normdGPS(idGPS);
        
        %% H
        HGalileo = [ndGalileo ones(size(ndGalileo,1),1)];
        HGPS = [ndGPS ones(size(ndGPS,1),1)];
        
        MGalileo = HGalileo' * HGalileo;%% 4x4
        MGPS = HGPS' * HGPS;%% 4x4
        
        HIntegration = [ndGalileo ones(size(ndGalileo,1),1) zeros(size(ndGalileo,1),1); ...
                        ndGPS zeros(size(ndGPS,1),1) ones(size(ndGPS,1),1)];%% since we have two independent clock errors.
        MIntegration = HIntegration'*HIntegration;
        %% PDOP
        PDOPGalileo = PDOP_calc(MGalileo);
        PDOPGPS = PDOP_calc(MGPS);
        PDOPIntegration = PDOP_calc(MIntegration);
        
        %% save
        times = [times;t];
        PDOPs = [PDOPs;[PDOPGalileo PDOPGPS PDOPIntegration]];
    end
    
    red_color = [153 0 0]/255;
    blue_color = [0 116 186]/255;
    orange_color = [223 80 35]/255;
    font_size = 16;
    fig = figure();
    set(fig,'defaulttextinterpreter','latex');
    plot(times,PDOPs(:,1),'LineStyle','--','LineWidth',2.5, 'Color', blue_color);grid on;hold on;
    plot(times,PDOPs(:,2),'LineStyle','-.','LineWidth',2.5, 'Color', orange_color);
    plot(times,PDOPs(:,3),'LineStyle','-','LineWidth',2.5, 'Color', red_color);
    title('PDOP over 24 hours','Interpreter','latex');
    xlabel('time: (s)','Interpreter','latex');
    ylabel('PDOP','Interpreter','latex');
    lgnd = legend({'Galileo Only','GPS Only','GPS+Galileo'}, 'Location', 'NorthWest');
    set(lgnd, 'Interpreter', 'Latex','FontSize', font_size);
    colormap summer
    set(gca,'TickLabelInterpreter','latex');
    xlabel('time: (s)','FontSize', font_size, 'Interpreter', 'latex');
    ylabel('PDOP','FontSize', font_size, 'Interpreter', 'latex');
    hold off;   
    
    
    
    
end

function [rec_llh, rec_xyz] = prepare()
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
end