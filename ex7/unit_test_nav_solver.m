function unit_test_nav_solver
%% source data from Allan Aasbjerg Nielsen
    close all;
    addpath('../utils')
    % true position (Landmaalervej, Hjortekaer)
    xtrue = [3507884.948 780492.718 5251780.403 0]';
    % positions of satellites 1, 4, 7, 13, 20, 24 and 25 in ECEF coordinate system, [m]
    xxyyzz = [16577402.072 5640460.750 20151933.185; ...
                11793840.229 -10611621.371 21372809.480; ...
                20141014.004 -17040472.264 2512131.115; ...
                22622494.101 -4288365.463 13137555.567; ...
                12867750.433 15820032.908 16952442.746; ...
                -3189257.131 -17447568.373 20051400.790; ...
                -7437756.358 13957664.984 21692377.935];
    pseudorange = [20432524.0 21434024.4 24556171.0 21315100.2 21255217.0 ...
                   24441547.2 23768678.3]'; % [m]
               
    sprior2 = 10^2; %5?2; %prior variance [m?2]
    
    options.usePrior = 0;options.useSOCP = 1;options.useDLT = 0;
    options.useWLS = 1;
    options.verbose = 0;
    options.x0_prior = [0 0 0 0]';
    options.maxiter = 100;
    options.threshold1 = 1e-6;
    options.threshold2 = 1e-6;
    options.prs_var = sprior2;
    res = navSolver(pseudorange, xxyyzz, options);
end