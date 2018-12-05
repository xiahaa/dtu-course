function ex_9_galileo_gps_display
    close all;
    clc;clear all;
    %% add utility functions
    addpath('../utils/')
    %% load Galileo Kepler Params
    GalileoKepParams;
    numGalileoSat = size(GalileoKeplerTb,1);    

    satGPS = {'A1','A2','A3','A4', ...
              'B1','B2','B3','B4', ...
              'C1','C2','C3','C4', ...
              'D1','D2','D3','D4', ...
              'E1','E2','E3','E4', ...
              'F1','F2','F3','F4'}';
          
    numGPSSat = size(satGPS,1);    

    %% orbit
    t = linspace(0,60*60*12,50);
    GalileosatPosX = zeros(numGalileoSat,numel(t));
    GalileosatPosY = zeros(numGalileoSat,numel(t));
    GalileosatPosZ = zeros(numGalileoSat,numel(t));

    for i = 1:numel(t)
        for j = 1:numGalileoSat
            [~, ~, ~, i1, i2, i3] = calc_sat_pos_with_Kepler(GalileoKeplerTb(j,:), 0, t(i));
            GalileosatPosX(j,i) = i1;
            GalileosatPosY(j,i) = i2;
            GalileosatPosZ(j,i) = i3;
        end
    end
    
    GPSsatPosX = zeros(numGPSSat,numel(t));
    GPSsatPosY = zeros(numGPSSat,numel(t));
    GPSsatPosZ = zeros(numGPSSat,numel(t));

    for i = 1:numel(t)
        for j = 1:numGPSSat
            slotID = satGPS{j};
            [~, ~, ~, i1, i2, i3] = calcSatPosition(slotID, 0, t(i));
            GPSsatPosX(j,i) = i1;
            GPSsatPosY(j,i) = i2;
            GPSsatPosZ(j,i) = i3;
        end
    end
    
    figure
    color = jet(9);
    markers = {'-o','*','s','d'};
    k = 1;
    for i = 1:numGalileoSat
        if mod(i,9) == 0
            h(k) = plot3(GalileosatPosX(i,:),GalileosatPosY(i,:),GalileosatPosZ(i,:),markers{1}, 'Color', color(k,:));hold on;
            k = k +1;
        else
            plot3(GalileosatPosX(i,:),GalileosatPosY(i,:),GalileosatPosZ(i,:),markers{1}, 'Color', color(k,:));hold on;
        end
    end
    
    for i = 1:numGPSSat
        if mod(i,4) == 0
            h(k) = plot3(GalileosatPosX(i,:),GalileosatPosY(i,:),GalileosatPosZ(i,:),markers{1}, 'Color', color(k,:));hold on;
            k = k + 1;
        else
            plot3(GPSsatPosX(i,:),GPSsatPosY(i,:),GPSsatPosZ(i,:),markers{2}, 'Color', color(k,:));hold on;
        end
    end
    [x,y,z] = sphere;
    x = x.*6371000;
    y = y.*6371000;
    z = z.*6371000;
    surf(x,y,z);
    grid on;
    title('Satellite Orbits','Interpreter','latex');
    legend([h(1), h(2), h(3), h(4), h(5), h(6), h(7), h(8), h(9)], ...
        'Galileo-1', 'Galileo-2', 'Galileo-3', 'GPS-A', 'GPS-B', ...
         'GPS-C', 'GPS-D',  'GPS-E', 'GPS-F');

    %% animation
    drawGif = 0;
    axis tight manual % this ensures that getframe() returns a consistent size
    giffilename = 'satellite.gif';
 
    t = linspace(0,60*60*12,12*60);
    GalileosatPosX = zeros(numGalileoSat,1);
    GalileosatPosY = zeros(numGalileoSat,1);
    GalileosatPosZ = zeros(numGalileoSat,1);
    GPSsatPosX = zeros(numGPSSat,1);
    GPSsatPosY = zeros(numGPSSat,1);
    GPSsatPosZ = zeros(numGPSSat,1);
    
    iter = 1;
    init = 0;
    frames = [];
    for i = 1:numel(t)
        for j = 1:numGalileoSat
            %% Galileo
            [~, ~, ~, i1, i2, i3] = calc_sat_pos_with_Kepler(GalileoKeplerTb(j,:), 0, t(i));
            GalileosatPosX(j,1) = i1;
            GalileosatPosY(j,1) = i2;
            GalileosatPosZ(j,1) = i3;
        end
        for j = 1:numGPSSat
            slotID = satGPS{j};
            [~, ~, ~, i1, i2, i3] = calcSatPosition(slotID, 0, t(i));
            GPSsatPosX(j,1) = i1;
            GPSsatPosY(j,1) = i2;
            GPSsatPosZ(j,1) = i3;
        end
        sats = [[GalileosatPosX GalileosatPosY GalileosatPosZ];[GPSsatPosX GPSsatPosY GPSsatPosZ]];
        if init == 0
            init = 1;
            satsh = ExampleHelperSat(numGalileoSat+numGPSSat,sats);
        else
            updatePlot(satsh, sats, t(i));
        end
        drawnow;
        frame = getframe(gcf); 
        frames = [frames;frame]; 
    end
    
    if drawGif == 1
        for i = 1:size(frames,1)
            im = frame2im(frames(i)); 
            [imind,cm] = rgb2ind(im,256); 
%             imwrite(imind,cm,strcat(num2str(i),'.png'));
            if iter == 1
                iter = 2;
                imwrite(imind,cm,giffilename,'gif', 'Loopcount',inf); 
            else
                imwrite(imind,cm,giffilename,'gif','WriteMode','append'); 
            end
        end
    end
    disp('Ending....');
end