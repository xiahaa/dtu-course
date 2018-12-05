function ex9_galileo_display
    close all;
    clc;clear all;

    addpath('../utils/')
    
    GalileoKepParams;
    numSat = size(GalileoKeplerTb,1);

    % Q1,2,3
    initTatPos = zeros(numSat,3);
    k = 1;
    for i = 1:numSat
        [~, ~, ~, i1, i2, i3] = calc_sat_pos_with_Kepler(GalileoKeplerTb(i,:), 0, 0);
        initTatPos(k,:) = [i1,i2,i3];
        k = k + 1;
    end
    figure
    color = jet(32);
    markers = ['o'];
    for i = 1:numSat
        plot3(initTatPos(i,1),initTatPos(i,2),initTatPos(i,3),'Marker',markers(1), ...
                'MarkerEdgeColor',color(i,:),'MarkerFaceColor',color(i,:),'MarkerSize', 8);hold on;
    end
    grid on;
    title('Satellite Positions','Interpreter','latex');
    
    %% Q4
    t = linspace(0,60*60*12,50);
    satPosX = zeros(numSat,numel(t));
    satPosY = zeros(numSat,numel(t));
    satPosZ = zeros(numSat,numel(t));

    for i = 1:numel(t)
        for j = 1:numSat
            [~, ~, ~, i1, i2, i3] = calc_sat_pos_with_Kepler(GalileoKeplerTb(j,:), 0, t(i));
            satPosX(j,i) = i1;
            satPosY(j,i) = i2;
            satPosZ(j,i) = i3;
        end
    end
    figure
    color = jet(32);
    markers = {'-o','*','s','d'};
    for i = 1:numSat
        plot3(satPosX(i,:),satPosY(i,:),satPosZ(i,:),markers{1}, 'Color', color(i,:));hold on;
    end
    [x,y,z] = sphere;
    x = x.*6371000;
    y = y.*6371000;
    z = z.*6371000;
    surf(x,y,z);
    grid on;
    title('Satellite Orbits','Interpreter','latex');
    
    drawGif = 0;
    axis tight manual % this ensures that getframe() returns a consistent size
    giffilename = 'satellite.gif';
 
    t = linspace(0,60*60*12,12*60);
    satPosX = zeros(numSat,1);
    satPosY = zeros(numSat,1);
    satPosZ = zeros(numSat,1);
    iter = 1;
    init = 0;
    frames = [];
    for i = 1:numel(t)
        for j = 1:numSat
            [~, ~, ~, i1, i2, i3] = calc_sat_pos_with_Kepler(GalileoKeplerTb(j,:), 0, t(i));
            satPosX(j,1) = i1;
            satPosY(j,1) = i2;
            satPosZ(j,1) = i3;
        end
        sats = [satPosX satPosY satPosZ];
        if init == 0
            init = 1;
            satsh = ExampleHelperSat(numSat,sats);
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
end






