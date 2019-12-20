function varargout = checkVisibility(numSat, dp, lat, lon)
    visibility = zeros(numSat,1);
    azimuths = zeros(numSat,1);
    zeniths = zeros(numSat,1);
    elevations = zeros(numSat,1);
    for i = 1:numSat
        %% conversion
        [e,n,u] = WGS842ENU(lat, lon, dp(i,1), dp(i,2), dp(i,3));
        %% compute azimuth and zenith
        [azi, zen, elevation] = calcAzimuthZenithElevation(e,n,u);
        if elevation > 5
            visibility(i) = 1;
        end
        azimuths(i) = azi;
        zeniths(i) = zen;
        elevations(i) = elevation;
    end
    varargout{1} = visibility;
    varargout{2} = azimuths;
    varargout{3} = zeniths;
    varargout{4} = elevations;
end