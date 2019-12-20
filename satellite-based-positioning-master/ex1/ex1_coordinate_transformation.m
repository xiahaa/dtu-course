function ex1_coordinate_transformation
    clear all;
    clc;
    consParams = struct('a',6378137.0,'f',1/298.257223563);
    
    %% add utils path
    addpath('../utils/')
    
    %% test 1, deci
    testCase(1,:) = [55.78575300466123,12.525384183973078,40];% DTU 101
    testCase(2,:) = [0,12.525384183973078,40];%equator
    testCase(3,:) = [90,0,40];
    testCase(4,:) = [-90,-15,40];
    
    testid = 1;% choose the desired point you want to test.
    
    lat = testCase(testid,1);% latitude
    lon = testCase(testid,2);% longitude
    height = testCase(testid,3);% altitude
    [x,y,z] = llhtoCartesian(lat, lon, height, consParams); % do forward transformation
    %% show
    disp(strcat('x: ',num2str(x),'; y: ',num2str(y),'; z: ',num2str(z)));
    %% reverse
    [latr, lonr, heightr] = Cartesian2llh(x,y,z,consParams);
    
    disp(strcat('true lat: ',num2str(lat),'; true lon: ',num2str(lon),'; true height: ',num2str(height)));% '; ellipsoid Height',num2str(heightr+N)));
    disp(strcat('lat: ',num2str(latr),'; lon: ',num2str(lonr),'; height: ',num2str(heightr)));% '; ellipsoid Height',num2str(heightr+N)));
    error = abs(lat-latr)+abs(lon-lonr)+abs(height-heightr);
    errorRMS = sqrt(mean(([lat, lon, height] - [latr, lonr, heightr]).^2));
    disp(strcat('error: ',num2str(error)));
    disp(strcat('RMS: ', num2str(errorRMS)));
    
    %% test case for deg minute second to decimal
    latd = 55;latm = 47;lats = 8.7108;
    lond = 12;lonm = 31;lons = 31.3824;
    lart = deg2deci(latd,latm,lats);
    lont = deg2deci(lond,lonm,lons);
    disp(strcat('lat: ', num2str(lat), '; latt', num2str(lart)));
    disp(strcat('lon: ', num2str(lon), '; lont', num2str(lont)));

end