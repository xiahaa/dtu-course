function ex4_interpolation_sat_position
    close all;
    clear all;
    
    % add utils
    addpath('../utils/');
    % sp3 file
    path = uigetfile('*.sp3','Open a SP3 file...');
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
    
    %% req 1: interpolation of 1 sat for 1 hour
    % given epoch
    epoch.year = 2018;
    epoch.Month = 9;
    epoch.Day = 17;
    
    %% find today's data
    dataOfToday = find_data(content, epoch);
    
    epoch.Hour = 12;
    epoch.Minute = 0;
    epoch.Second = 0;
    satId = 5;
    
    keypts = [];
    cubic_sat_trace = [];
    lagrange_sat_trace = [];
    for i = 0:1:60
        % every minute
        epoch.Minute = epoch.Minute + 1;
        %% find the nearest 5 data in today's data
        neighbor_ids = find_neighbor_ids(dataOfToday, epoch, 2);
        %% interpolation 1 time
        t_interp = epoch.Hour*3600+epoch.Minute*60+ epoch.Second;
        [sx,sy,sz] = interp_sat_pos(satId, neighbor_ids, dataOfToday, t_interp, '1');
        cubic_sat_trace = [cubic_sat_trace;[sx sy sz]];
        [sx,sy,sz] = interp_sat_pos(satId, neighbor_ids, dataOfToday, t_interp, '2');
        lagrange_sat_trace = [lagrange_sat_trace;[sx sy sz]];
        if mod(i,15) == 0
            keypts = [keypts;[sx sy sz]];
        end
        
    end
    
    figure(1)
    h1 = plot3(cubic_sat_trace(:,1),cubic_sat_trace(:,2),cubic_sat_trace(:,3),'r-','LineWidth',2);
    grid on; hold on;
    h3 = plot3(lagrange_sat_trace(:,1),lagrange_sat_trace(:,2),lagrange_sat_trace(:,3),'m-','LineWidth',2);
    h2 = plot3(keypts(:,1),keypts(:,2),keypts(:,3),'bs','MarkerSize',8);
    legend([h1,h2,h3],'1 hour cubic spline interpolation: trace', 'Precise Position from SP3', ...
        '1 hour lagrange interpolation: trace', 'Interpreter','latex');
    xlabel('x:(m)','Interpreter','latex');ylabel('y:(m)','Interpreter','latex');zlabel('z:(m)','Interpreter','latex');
    title(strcat('Position interpolation for sat ', num2str(satId)),'Interpreter','latex');
    
    %% req 2: remove 1, interpolate and compare
    % given epoch
    epoch2.year = 2018;
    epoch2.Month = 9;
    epoch2.Day = 17;
    epoch2.Hour = 12;
    epoch2.Minute = 0;
    epoch2.Second = 0;
    satId = 5;
    
    %% find today's data
    dataOfToday2 = find_data(content, epoch2);    
    cnt = 2;
    neighbor_ids = find_neighbor_ids(dataOfToday2, epoch2, cnt);
    
    err_cubic_spline = zeros(2*cnt+1,1);
    for i = 1:numel(err_cubic_spline)
        t_interp = dataOfToday2{neighbor_ids(i)}.Hour*3600+dataOfToday2{neighbor_ids(i)}.Minute*60+ dataOfToday2{neighbor_ids(i)}.Second;
        newid = neighbor_ids;
        newid(i) = [];
        tsx = dataOfToday2{neighbor_ids(i)}.satPos(satId).x;
        tsy = dataOfToday2{neighbor_ids(i)}.satPos(satId).y;
        tsz = dataOfToday2{neighbor_ids(i)}.satPos(satId).z;
        [sx,sy,sz] = interp_sat_pos(satId, newid, dataOfToday2, t_interp, '1');
        err_cubic_spline(i) = norm([tsx-sx,tsy-sy,tsz-sz]);
    end
    figure(2);
    colormap(cool);
    bar([err_cubic_spline],'grouped');
    hold on;grid on;
    xlabel('x: interpolation id', 'Interpreter','latex');
    ylabel('y: interpolation error (m)', 'Interpreter','latex');
    title('Error Analysis','Interpreter','latex');
    
    %% req 3, a list of satellite position at 12:05
    % given epoch
    epoch3.year = 2018;
    epoch3.Month = 9;
    epoch3.Day = 17;
    epoch3.Hour = 12;
    epoch3.Minute = 5;
    epoch3.Second = 0;
    
    cnt = 2;
    
    epoch4.year = 2018;
    epoch4.Month = 9;
    epoch4.Day = 17;
    epoch4.Hour = 12;
    epoch4.Minute = 0;
    epoch4.Second = 0;
    
    satPosPred = zeros(content.satNum,3);
    satPosSP3 = zeros(content.satNum,3);
    t_interp = epoch3.Hour*3600+epoch3.Minute*60+ epoch3.Second;
    dataOfToday3 = find_data(content, epoch3);    
    neighbor_ids = find_neighbor_ids(dataOfToday3, epoch3, cnt);
    closet_id = find_closet_id(dataOfToday3,epoch4);
    for i = 1:content.satNum
        satId = i;
        [sx,sy,sz] = interp_sat_pos(satId, neighbor_ids, dataOfToday3, t_interp, '1');
        satPosPred(i,:) = [sx,sy,sz];
        tsx = dataOfToday3{closet_id}.satPos(satId).x;
        tsy = dataOfToday3{closet_id}.satPos(satId).y;
        tsz = dataOfToday3{closet_id}.satPos(satId).z;
        satPosSP3(i,:) = [tsx,tsy,tsz];
    end
    disp(satPosPred);
    figure(3);
    [xe,ye,ze] = sphere;
    xe = xe.*6371000;
    ye = ye.*6371000;
    ze = ze.*6371000;
    surf(xe,ye,ze);
    grid on;
    hold on;
    h4=plot3(satPosSP3(:,1),satPosSP3(:,2),satPosSP3(:,3),'rs','MarkerSize',5);
    h5=plot3(satPosPred(:,1),satPosPred(:,2),satPosPred(:,3),'bd','MarkerSize',5);
    for i = 1:content.satNum
        text(satPosPred(i,1),satPosPred(i,2),satPosPred(i,3)+2000,num2str(i),'FontSize',15);
    end
    title('Satellite Position','Interpreter','latex');
    legend([h4,h5],'Positions at 12:00', 'Positions at 12:05', 'Interpreter','latex');
    xlabel('x:(m)','Interpreter','latex');ylabel('y:(m)','Interpreter','latex');zlabel('z:(m)','Interpreter','latex');
    
    %% req 4, estimate movement
    % given epoch
    epoch5.year = 2018;
    epoch5.Month = 9;
    epoch5.Day = 17;
    epoch5.Hour = 12;
    epoch5.Minute = 0;
    epoch5.Second = 0;
    
    cnt = 2;
    sampleNum = 500;
    ts = linspace(0,5*60,sampleNum)';
    satPosPred2 = zeros(content.satNum, numel(ts), 3);
    dataOfToday5 = find_data(content, epoch5);    
    for i = 1:numel(ts)
        t_interp = epoch5.Hour*3600+ts(i);
        neighbor_ids = find_neighbor_ids_with_ts(dataOfToday5, t_interp, cnt);
        for j = 1:content.satNum
            satId = j;
            [sx,sy,sz] = interp_sat_pos(satId, neighbor_ids, dataOfToday5, t_interp, '1');
            satPosPred2(j, i, :) = [sx,sy,sz];
        end
    end
    %% numerical integration
    moves = zeros(content.satNum, 1);
    for j = 1:content.satNum
        satPosJ = satPosPred2(j,:,:);
        satPosJ = reshape(satPosJ,sampleNum,3,1);
        posdiff = satPosJ(2:end,:) - satPosJ(1:end-1,:);
        moves(j) = sum(sqrt(posdiff(:,1).^2+posdiff(:,2).^2+posdiff(:,3).^2));
    end
    figure(4);
    colormap(cool);
    bar([moves],'grouped');
    hold on;grid on;
    xlabel('x: satellite id', 'Interpreter','latex');
    ylabel('y: movement (m)', 'Interpreter','latex');
    title('Movements in 5 minutes','Interpreter','latex');
    
    slotID = strcat('A',num2str(1));
    t0 = 12*3600;
    t = 12*3600+60*5;
    s = move_estimation_via_RK(slotID, t0, t);
    disp(s);
end






