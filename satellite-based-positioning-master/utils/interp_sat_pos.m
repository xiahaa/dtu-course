function [sx,sy,sz] = interp_sat_pos(satId, neighbor_ids, dataOfToday, t_interp, interp_opt)
    % neighboring data
    format long;
    satpos = zeros(numel(neighbor_ids),4);
    for i = 1:1:numel(neighbor_ids)
        t = dataOfToday{neighbor_ids(i)}.Hour*3600+dataOfToday{neighbor_ids(i)}.Minute*60+ ...
            dataOfToday{neighbor_ids(i)}.Second;
        p = [dataOfToday{neighbor_ids(i)}.satPos(satId).x, ...
             dataOfToday{neighbor_ids(i)}.satPos(satId).y, ...
             dataOfToday{neighbor_ids(i)}.satPos(satId).z];
        satpos(i,:) = [t, p];
    end
    
    if interp_opt == '1'
        % interpolation
        sx = interpolation1D(satpos(:,1),satpos(:,2),t_interp,'spline');
        sy = interpolation1D(satpos(:,1),satpos(:,3),t_interp,'spline');
        sz = interpolation1D(satpos(:,1),satpos(:,4),t_interp,'spline');
    elseif interp_opt == '2'
        % interpolation
        sx = lagrange_interpolation(satpos(:,1),satpos(:,2),t_interp);
        sy = lagrange_interpolation(satpos(:,1),satpos(:,3),t_interp);
        sz = lagrange_interpolation(satpos(:,1),satpos(:,4),t_interp);
    end
end

function y = interpolation1D(sx,sy,x,method)
    y = interp1(sx,sy,x,method);
end