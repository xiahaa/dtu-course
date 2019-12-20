function [prs_full, sat_pos_full, visibilities, xsol] = simulation7(epoch, content, llh, TECU, recerr)
    %% add utils path
    addpath('../utils/');
    addpath('../utils/3rdparty');
    
    x0 = [0;0;0];
    x0 = [x0;0];%% add initial receiver clock error
    
    %% options for the navSolver
    %   for initialization:
    %       1: use the x0 provided by the user;
    %       2: use the proposed second-order-cone programming;
    %       3: use the proposed direct linear transformation method;
    %       4: use the Bancraft method;
    %   for iterative solver:
    %       1: Gauss-Newton;
    %       2: Steepest Dscent;
    %       3: Levenberg-Marquardt;   
    options.initialization = 1;% can be 1,2,3,4
    options.solver = 1;% can be 1,2,3
    options.verbose = 0;% will generate intermediate print info
    options.maxiter = 100;% maximum iteration number
    %% two thresholds for iteration: smaller, more iteration.
    options.threshold1 = 1e-6;%  
    options.threshold2 = 1e-6;%
    
    %% generate test case
    [prs, sat_pos, visibilities, sat_pos_full, prs_full] = find_sat_xyz(epoch,llh, content, TECU, recerr);
    % pseudoranges with atomospheric delays and satellite clock error
    % corrected
    pseudorange_full_cor = correct_error_excep_recever(prs(:,1), prs(:,3), prs(:,4), prs(:,5));

    options.x0_prior = x0;
    %% solve with fully corrected pseudorange
    sprior2 = 10^2; %5^2; %prior variance [m^2]
    options.prs_var = sprior2;% prior covariance
    [x_cor,std_x_cor,QDOP_cor,Qenu_cor,llh3] = navSolver(pseudorange_full_cor, sat_pos, options);
    PDOP_cor = sqrt(trace(QDOP_cor(1:3,1:3)));
    GDOP_cor = sqrt(trace(QDOP_cor(1:4,1:4)));
    HDOP_cor = sqrt(trace(Qenu_cor(1:2,1:2)));
    VDOP_cor = sqrt(trace(Qenu_cor(3,3)));
    
    xsol = x_cor;
end

function prs = correct_error_excep_recever(prs_raw, d_iono, d_trop, d_satclk)
%% correct_error_excep_recever: 
%   correct pseudoranges with estimated atmospeheric effects.
%   inputs: 
%       prs_raw: raw pseudoranges;
%       d_iono: ionospheric delay;
%       d_trop: tropospheric delay;
%       d_satclk: satellite clock error;
%   outputs:
%       prs: corrected pseudoranges
%% Author: xiahaa@space.dtu.dk

    prs = prs_raw;
    if ~isempty(d_iono)
        prs = prs - d_iono;
    end
    if ~isempty(d_trop)
        prs = prs - d_trop;
    end
    if ~isempty(d_satclk)
        prs = prs + d_satclk;
    end
end

function [prs, sat_pos, visibilities, sat_pos_full, prs_full] = find_sat_xyz(epoch,rec_llh, content, TECU, recerr)
%% generate_test_case1: 
%   generate pseudoranges for a hard-coded postision (DTU 101).
%   inputs: none
%   outputs:
%       prs: pseudoranges
%       sat_pos: satellite positions
%       x0: true receiver position
%% Author: xiahaa@space.dtu.dk

    %% read tpday's information
    dataOfToday = find_data(content, epoch);   
    
    %% satellite position prediction
    satPosPred = zeros(content.satNum,3);
    t_interp = epoch.Hour*3600+epoch.Minute*60+ epoch.Second;
    cnt = 2;
    neighbor_ids = find_neighbor_ids(dataOfToday, epoch, cnt);
    for i = 1:content.satNum
        satId = i;
        [sx,sy,sz] = interp_sat_pos(satId, neighbor_ids, dataOfToday, t_interp, '1');
        satPosPred(i,:) = [sx,sy,sz];
    end
    
    %% parse content
    satInfo = cell(content.satNum,1);
    for i = 1:size(dataOfToday,1)
        data = dataOfToday{i,1};
        time = data.Hour * 3600 + data.Minute * 60 + data.Second;
        for j = 1:size(data.satPos,1)
            satinfo = [time, data.satPos(j).x, data.satPos(j).y, data.satPos(j).z,data.satPos(j).clock];
            satInfo{j}=[satInfo{j};satinfo];
        end
    end 
    
    %%
    ts = t_interp;
    prs = [];
    ids = [];
    visibilities = zeros(content.satNum,1);
    prs_full = zeros(content.satNum,1);
    for i = 1:content.satNum
        [visibility, pr, R, d_iono, d_trop, d_satclk, d_recclk, clkerr] = calc_pseudorange(rec_llh(1,1), rec_llh(1,2), rec_llh(1,3), ...
                                      satPosPred(i,:), satInfo, i, ts, TECU, recerr);
        if visibility == 1
            visibilities(i) = 1;
            ids = [ids;i];
            prs = [prs;[pr, R, d_iono, d_trop, d_satclk, d_recclk, clkerr]];
            prs_full(i) = pr;
        end
    end
       
    %% assign output
    sat_pos = satPosPred(ids,:);
    sat_pos_full = satPosPred;
    
end
