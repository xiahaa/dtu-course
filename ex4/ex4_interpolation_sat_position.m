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
    % given epoch
    epoch.year = 2018;
    epoch.Month = 9;
    epoch.Day = 17;
    epoch.Hour = 12;
    epoch.Minute = 5;
    epoch.Second = 0;
    satId = 5;
    % find today's data
    dataOfToday = {};
    j = 1;
    for i = 1:size(content.sections,1)
        timediff = abs(content.sections{i}.year - epoch.year) * 365 + ...
                   abs(content.sections{i}.Month - epoch.Month) * 30 + ...
                   abs(content.sections{i}.Day - epoch.Day);
        % not today, continue
        if timediff >= 1
            continue;
        end
        % today's data
        da.Hour = content.sections{i}.Hour;
        da.Minute = content.sections{i}.Minute;
        da.Second = content.sections{i}.Second;
        da.satPos = content.sections{i}.satPos;
        dataOfToday{j,1} = da;
        j = j + 1;
    end
    
    % find the nearest 5 data in today's data
    clockErr = [];
    for i = 1:size(dataOfToday,1)
        clockErr(i) = (dataOfToday{i}.Hour - epoch.Hour)*3600 + ...
                      (dataOfToday{i}.Minute - epoch.Minute)*60 + ...
                      (dataOfToday{i}.Second - epoch.Second);
    end    
    [minval, minid] = min(abs(clockErr));
    neightborIDs = [];
    if (minid - 2) >= 1 && (minid + 2) <= numel(clockErr)
        neightborIDs = (minid - 2):1:(minid + 2);
    else
        if (minid - 2) >= 1
            neightborIDs = (minid - 4):1:minid;
        else
            neightborIDs = minid:1:(minid+4);
        end
    end
    % neighboring data
    satpos = zeros(numel(neightborIDs),4);
    for i = 1:1:numel(neightborIDs)
        t = dataOfToday{neightborIDs(i)}.Hour*3600+dataOfToday{neightborIDs(i)}.Minute*60+ ...
            dataOfToday{neightborIDs(i)}.Second;
        p = [dataOfToday{neightborIDs(i)}.satPos(satId).x, ...
             dataOfToday{neightborIDs(i)}.satPos(satId).y, ...
             dataOfToday{neightborIDs(i)}.satPos(satId).z];
        satpos(i,:) = [t, p];
    end
    % interpolation
    intert = epoch.Hour*3600+epoch.Minute*60+ epoch.Second;
    sx = interpolation1D(satpos(:,1),satpos(:,2),intert,'spline');
    sy = interpolation1D(satpos(:,1),satpos(:,3),intert,'spline');
    sz = interpolation1D(satpos(:,1),satpos(:,4),intert,'spline');
    
    figure(1)
    ts = linspace(satpos(1,1),satpos(end,1),1000);
    spilinx = interpolation1D(satpos(:,1),satpos(:,2),ts,'spline');
    spiliny = interpolation1D(satpos(:,1),satpos(:,3),ts,'spline');
    spilinz = interpolation1D(satpos(:,1),satpos(:,4),ts,'spline');
    plot3(spilinx,spiliny,spilinz,'g-','LineWidth',2);
    hold on; grid on;
    plot3(satpos(:,2),satpos(:,3),satpos(:,4),'rd','MarkerSize',5);
    plot3(sx,sy,sz,'bs','MarkerSize',5);
    title("Interpolation of Satellite's positions",'Interpreter','latex');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    legend({'Spline','SP3 Sat Pos', 'Interpolated Sat Pos'},'Interpreter','latex');
end

function y = interpolation1D(sx,sy,x,method)
    y = interp1(sx,sy,x,method);
end


