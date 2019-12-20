function ex3_enu_conversion

    addpath('../utils/');

    % fetched satellite position from igs at *  2018  9 18  5 30  0.00000000
    data = [-14084.965421 -13335.977939 -18444.389348; ... 
            -8412.178797  14931.928610  20810.047600; ...
            -13642.570569 -22690.514868   2505.504310; ...
            -390.404350  21748.781259  14982.038597; ...
            -22918.304145   6890.582919  11571.183052; ...
            -25196.382195  -6397.745158   6768.549040; ...
            2796.631366 -22395.316629 -13830.185289; ...
            -13947.715293  -6618.056480  21577.959060; ...
            15947.740945  -1338.501133 -21183.204680; ...
            15947.740945  -1338.501133 -21183.204680; ...
            9161.084379  24938.647725  -2378.996417; ...
            -11093.440122  22608.322634  -8521.726915; ...
            18485.647522 -18026.093616  -5511.098696; ...
            3157.925831  20557.837992 -16521.290198; ...
            3402.684664 -20247.785807  16449.616936; ...
            -19351.327692  10977.974490 -14018.876526; ...
            -4035.301059 -15495.645056 -21685.320944; ...
            -19015.304410  17811.742561  -5797.886412; ...
            19708.841479   7637.705325 -16267.555982; ...
            26132.764030   2961.941803   4275.236351; ...
            -9918.277161 -24070.081964  -5110.381059; ...
            -7975.048024 -15242.870677  20435.202393; ...
            9811.195637  14662.117690 -20056.505271; ...
            16842.930627  18008.705148   9990.949356; ...
            10105.122386 -12626.042864  21040.684542; ...
            11417.769712 -23432.342744  -4392.131096; ...
            -15790.334354  -2544.884817 -20573.112489; ...
            12353.007291   7962.709748  22097.350060; ...
            -26312.156995    524.630853  -4199.834644; ...
            21865.789063  -7671.878821  13389.204161; ...
            19113.591091 -12285.821643 -13637.262853; ...
            ];
    data = data.*1000;% km to m
    % local position llh
    testCase(1,:) = [55.78575300466123,12.525384183973078,0];% DTU 101
    consParams = struct('a',6378137.0,'f',1/298.257223563); % some constants
    testid = 1;
    lat = testCase(testid,1);%
    lon = testCase(testid,2);%
    height = testCase(testid,3);%
    [xo,yo,zo] = llhtoCartesian(lat, lon, height, consParams);% to ECEF
    
    % satpos in enu
    satPosENU = zeros(size(data,1),3);
    azimuths = zeros(size(data,1),1);
    zeniths = zeros(size(data,1),1);
    visibles = zeros(size(data,1),1);
    dists = ones(size(data,1),1).*-1;
    for i = 1:size(data,1)
        satPos = data(i,:);
        dx = satPos(1) - xo;
        dy = satPos(2) - yo;
        dz = satPos(3) - zo;
        % conversion
        [e,n,u] = WGS842ENU(lat, lon, dx, dy, dz);
        satPosENU(i,:) = [e,n,u];
        
        %% compute azimuth and zenith
        [azimuth, zenith, elevation] = calcAzimuthZenithElevation(e,n,u);
        azimuths(i) = azimuth;
        zeniths(i) = zenith;

        if elevation > 5 % visibility threshold
            % identify as visible
            visibles(i) = 1;
            dists(i) = sqrt(dx^2+dy^2+dz^2);
        else
            visibles(i) = 0;
        end
    end
    
    % earth
    [xe,ye,ze] = sphere;
    xe = xe.*6371000;
    ye = ye.*6371000;
    ze = ze.*6371000;
    surf(xe,ye,ze);
    grid on;
    hold on;
    
    % local point
    plot3(xo,yo,zo,'m*','MarkerSize',10);
    
    for i = 1:size(data,1)
        if visibles(i) == 1
            % shown as a diamond with red color
            plot3(data(i,1),data(i,2),data(i,3),'rd','MarkerSize',10);
            % line satellite to local point
            line = [[xo, yo, zo];[data(i,1),data(i,2),data(i,3)]];
            plot3(line(:,1),line(:,2),line(:,3),'LineWidth',3);
            % show the distance
            text(0.5*(xo+data(i,1)),0.5*(yo+data(i,2)),0.5*(zo+data(i,3)),strcat('dist=',num2str(dists(i))),'Interpreter','latex');
        else
            % if not visible, display as black square.
            plot3(data(i,1),data(i,2),data(i,3),'ks','MarkerSize',10);
        end
    end
    title('Assignment3','Interpreter','latex');
end