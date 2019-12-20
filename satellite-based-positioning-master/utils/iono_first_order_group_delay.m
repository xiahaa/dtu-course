function d_iono = iono_first_order_group_delay(ts, tecfile, el)

    if isempty(tecfile) == 1
        tecfile = '../data/tec_2018_09_16.mat';
    end
    %% load relevent tec file;
    % file can be downloaded from: http://sesolstorm.kartverket.no/moreplots.xhtml
    teclut = readteclut(tecfile);

    %% find corresponding tec
    timediff = abs(teclut(:,1) - ts);
    [minval, minid] = min(timediff);
    interpids = [];
    if minid == 1
        interpids = [1 2 3]';
    elseif minid == size(teclut,1)
        interpids = [minid - 2, minid - 1, minid]';
    else
        interpids = [minid-1, minid, minid + 1]';
    end
    %% fetch tec values from southern Norway
    mes = teclut(interpids,2);
    tecm=interp1(interpids,mes,ts,'cubic');
    
    f1 = 1.57542e9;%% L1 frequency
    RE = 6371e3;
    hI = 350e3;
    d_ionov = 40.3/f1^2*tecm;%%omit 1e16 since tecm is normalized
    
    %% obliquity factor
    zenith = pi/2 - el;
    OF = (1-((RE*sin(zenith))/(RE+hI))^2)^(-1/2);
    d_iono = d_ionov * OF;
end

function teclut = readteclut(tecfile)
%% read TECU data
    teclut=load(tecfile); 
end