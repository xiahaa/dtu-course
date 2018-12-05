function hpol = skyplot_lite(azi1, el1, prn1, color1, marker1, ...
                          azi2, el2, prn2, color2, marker2)
    %% Prepare axis =========================================================== 
    % hAxis = newplot(hAxis); 
    hAxis = gca;

    %--- Get x-axis text color so grid is in same color ----------------------- 
    tc = get(hAxis, 'xcolor'); 

    hold(hAxis, 'on'); 

    %--- Plot white background ------------------------------------------------ 
    rectangle('position', [-90, -90, 180, 180], ... 
              'Curvature', [1 1], ... 
              'facecolor', 'white', ... 
              'edgecolor', tc); 
 
    %% Plot spokes ============================================================ 

    %--- Find spoke angles ---------------------------------------------------- 
    % Only 6 lines are needed to divide circle into 12 parts 
    th = (1:6) * 2*pi / 12; 

    %--- Convert spoke end point coordinate to Cartesian system --------------- 
    cst = cos(th); snt = sin(th); 
    cs = [cst; -cst]; 
    sn = [snt; -snt]; 

    %--- Plot the spoke lines ------------------------------------------------- 
    line(90*sn, 90*cs, 'linestyle', ':', 'color', tc, 'linewidth', 0.5, ... 
        'handlevisibility', 'off'); 

    %% Annotate spokes in degrees ============================================= 
    rt = 1.1 * 90; 

    for i = 1:max(size(th)) 

        %--- Write text in the first half of the plot ------------------------- 
        text(rt*snt(i), rt*cst(i), int2str(i*30), ... 
            'horizontalalignment', 'center', 'handlevisibility', 'off'); 

        if i == max(size(th)) 
            loc = int2str(0); 
        else 
            loc = int2str(180 + i*30); 
        end 

        %--- Write text in the opposite half of the plot ---------------------- 
        text(-rt*snt(i), -rt*cst(i), loc, ... 
            'handlevisibility', 'off', 'horizontalalignment', 'center'); 
    end 

    %% Plot elevation grid ==================================================== 

    %--- Define a "unit" radius circle ---------------------------------------- 
    th = 0 : pi/50 : 2*pi; 
    xunit = cos(th); 
    yunit = sin(th); 

    %--- Plot elevation grid lines and tick text ------------------------------ 
    for elevation = 0 : 30 : 90 
        elevationSpherical = 90*cos((pi/180) * elevation); 

        line(yunit * elevationSpherical, xunit * elevationSpherical, ... 
            'lineStyle', ':', 'color', tc, 'linewidth', 0.5, ... 
            'handlevisibility', 'off'); 

        text(0, elevationSpherical, num2str(elevation), ... 
            'BackgroundColor', 'white', 'horizontalalignment','center', ... 
            'handlevisibility', 'off'); 
    end 

    %--- Set view to 2-D ------------------------------------------------------ 
    view(0, 90); 

    %--- Set axis limits ------------------------------------------------------ 
    %save some space for the title 
    axis([-95 95 -90 101]); 

    %% Transform elevation angle to a distance to the center of the plot ------ 
    elSpherical1 = 90*cos(el1 * pi/180); 
    elSpherical2 = 90*cos(el2 * pi/180); 
    %--- Transform data to Cartesian coordinates ------------------------------ 
    yy1 = elSpherical1 .* cos(azi1 * pi/180); 
    xx1 = elSpherical1 .* sin(azi1 * pi/180); 
    yy2 = elSpherical2 .* cos(azi2 * pi/180); 
    xx2 = elSpherical2 .* sin(azi2 * pi/180); 
    
    %% Plot data on top of the grid =========================================== 
    hpol = plot(hAxis, xx1', yy1', marker1, 'Color', color1, 'MarkerSize', 10, ...
        'MarkerEdgeColor', color1, 'MarkerFaceColor', color1); hold on;
    hpol = plot(hAxis, xx2', yy2', marker2, 'Color', color2, 'MarkerSize', 10, ...
        'MarkerEdgeColor', color2, 'MarkerFaceColor', color2); hold on;
    
    %--- Place satellite PRN numbers at the latest position ------------------- 
    for i = 1:length(prn1) 
        if(prn1(i) ~= 0) 
            % The empthy space is used to place the text a side of the last 
            % point. This solution results in constant offset even if a zoom 
            % is used. 
            text(xx1(i, end), yy1(i, end), ['  ', int2str(prn1(i))], 'color', color1); 
        end 
    end 

    for i = 1:length(prn2) 
        if(prn2(i) ~= 0) 
            % The empthy space is used to place the text a side of the last 
            % point. This solution results in constant offset even if a zoom 
            % is used. 
            text(xx2(i, end), yy2(i, end), ['  ', int2str(prn1(i))], 'color', color2); 
        end 
    end 
    
    %--- Make sure both axis have the same data aspect ratio ------------------ 
    axis(hAxis, 'equal'); 

    %--- Switch off the standard Cartesian axis ------------------------------- 
    axis(hAxis, 'off'); 
end