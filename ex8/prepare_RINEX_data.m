function varargout = prepare_RINEX_data(contentRINEX, date)
%% prepare_RINEX_data: extract pseudoranges and satellite id
% Author: xiahaa@space.dtu.dk
    data = extract_data(contentRINEX, date);
    [prs, satIDs] = extract_pr(data);
    varargout{1} = data;
    varargout{2} = prs;
    varargout{3} = satIDs;
end

function varargout = extract_pr(data)
%% here, extract the pseudoranges (L1-L2) from RINEX file.
%  todo: whether or not we need to correct for the Doppler shift.
% Author: xiahaa@space.dtu.dk
    prs = zeros(size(data.SNRs,1),1);
    satIDs = zeros(size(data.SNRs,1),1);
    for i = 1:size(data.SNRs,1)
        prs(i) = data.satPos(i).C1;
        satIDs(i) = str2num(data.SNRs{i}(2:end));
    end
    varargout{1} = prs;
    varargout{2} = satIDs;
end

function varargout = extract_data(content, date)
%% extract data corresponding to the given date from the raw and complete
%  RINEX content.
% Author: xiahaa@space.dtu.dk
    dy = content.date.year - date.year;
    dm = content.date.Month - date.Month;
    dd = content.date.Day - date.Day;
    if abs(dy)>1 || abs(dm)>1 || abs(dd)>1
        error('The date does not correspond to the given RINEX file!!!!');
        varargout{1} = [];
        return;
    end
    ts1 = date.Hour * 3600 + date.Minute * 60 + date.Second;
    nearest = -1;
    nesrestDt = 1e6;
    for i = 1:size(content.sections,1)
        section = content.sections{i};
        ts2 = section.Hour * 3600 + section.Minute * 60 + section.Second;

        if abs(ts2 - ts1) < nesrestDt
            nesrestDt = abs(ts2 - ts1);
            nearest = i;
        end
    end
    
    if nesrestDt < 1
        data = content.sections{nearest};
        varargout{1} = data;
        disp('found!');
    else
        data = content.sections{nearest};
        varargout{1} = data;
        warning('Only found the nearest!');
    end
end