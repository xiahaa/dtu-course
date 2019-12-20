function varargout = prepare_SP3_data(contentSP3, epoch)
%% prepare_SP3_data: extract satellite position and satellite id
% Author: xiahaa@space.dtu.dk
    addpath('../utils/');
    dataOfToday = find_data(contentSP3, epoch);
    satInfo = parse_sat_info(dataOfToday, contentSP3.satNum);
    sp3SatIDs = extract_sat_IDs(contentSP3.satNames, contentSP3.satNum);
    varargout{1} = dataOfToday;
    varargout{2} = satInfo;
    varargout{3} = sp3SatIDs;
end

function satInfo = parse_sat_info(dataOfToday, satNum)
%% extract satellite position and clock error
% Author: xiahaa@space.dtu.dk
    satInfo = cell(satNum,1);
    for i = 1:size(dataOfToday,1)
        data = dataOfToday{i,1};
        time = data.Hour * 3600 + data.Minute * 60 + data.Second;
        for j = 1:size(data.satPos,1)
            satinfo = [time, data.satPos(j).x, data.satPos(j).y, data.satPos(j).z,data.satPos(j).clock];
            satInfo{j}=[satInfo{j};satinfo];
        end
    end
end

function sp3SatIDs = extract_sat_IDs(sat_names, sat_num)
%% extract satellite id in sp3 file.
% Author: xiahaa@space.dtu.dk
    for i = 1:sat_num
        sp3SatIDs(i) = str2num(sat_names(i,2:end));
    end
end