function ex2_satellite_position_calc
    close all;
    clc;clear all;
    orbitName = ['A','B','C','D','E','F'];% satellite orbital name
    satId = [1,2,3,4];% satellite id
    
    addpath('../utils/')
    
    % Q1,2,3
    initTatPos = zeros(numel(orbitName)*numel(satId),3);
    k = 1;
    for i = 1:numel(orbitName)
        for j = 1:numel(satId)
            slotID = strcat(orbitName(i),num2str(satId(j)));% satellite name
            [q1, q2, q3, i1, i2, i3] = calcSatPosition(slotID, 0, 0);% positions for given satellite
            initTatPos(k,:) = [i1,i2,i3];
            k = k + 1;
        end
    end
    figure
    color = jet(24);
    markers = ['o','*','s','d'];
    for i = 1:numel(orbitName)
        for j = 1:numel(satId)
            k = (i - 1)*numel(satId)+j;
            plot3(initTatPos(k,1),initTatPos(k,2),initTatPos(k,3),'Marker',markers(j),'Color',color(i,:));hold on;
        end
    end
    grid on;
    title('Satellite Positions','Interpreter','latex');
    
    %% Q4
    t = linspace(0,60*60*12,50);
    satPosX = zeros(numel(orbitName)*numel(satId),numel(t));
    satPosY = zeros(numel(orbitName)*numel(satId),numel(t));
    satPosZ = zeros(numel(orbitName)*numel(satId),numel(t));

    for i = 1:numel(t)
        for j = 1:numel(orbitName)
            for k = 1:numel(satId)
                slotID = strcat(orbitName(j),num2str(satId(k)));
                [q1, q2, q3, i1, i2, i3] = calcSatPosition(slotID, 0, t(i));
                satPosX((j-1)*numel(satId)+k,i) = i1;
                satPosY((j-1)*numel(satId)+k,i) = i2;
                satPosZ((j-1)*numel(satId)+k,i) = i3;
            end
        end
    end
    figure
    color = jet(24);
    markers = ['o','*','s','d'];
    for i = 1:numel(orbitName)
        for j = 1:numel(satId)
            k = (i - 1)*numel(satId)+j;
            plot3(satPosX(k,:),satPosY(k,:),satPosZ(k,:),'Marker',markers(j),'Color',color(k,:));hold on;
        end
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
    satPosX = zeros(numel(orbitName)*numel(satId),1);
    satPosY = zeros(numel(orbitName)*numel(satId),1);
    satPosZ = zeros(numel(orbitName)*numel(satId),1);
    iter = 1;
    init = 0;
    frames = [];
    for i = 1:numel(t)
        for j = 1:numel(orbitName)
            for k = 1:numel(satId)
                slotID = strcat(orbitName(j),num2str(satId(k)));
                [q1, q2, q3, i1, i2, i3] = calcSatPosition(slotID, 0, t(i));
                satPosX((j-1)*numel(satId)+k,1) = i1;
                satPosY((j-1)*numel(satId)+k,1) = i2;
                satPosZ((j-1)*numel(satId)+k,1) = i3;
            end
        end
        sats = [satPosX satPosY satPosZ];
        if init == 0
            init = 1;
            satsh = ExampleHelperSat(24,sats);
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






