function varargout = RINEX_file_parser(path)
%% parsing a RINEX file and extract useful information.
% Author: xiahaa@space.dtu.dk
    %% fopen file
    fid = fopen(path,"r");
    linenum = 1;
    satNum = 0;
    shiftbase = 25;

    sections = {};
    k = 1;
    kk = 1;
    satPos =struct('L1',0,'L2',0,'C1',0,'P2',0,'P1',0,'D1',0,'D2',0);
    
    while ~feof(fid)
        tline = fgetl(fid);
        if linenum == 13
            %% extract the approx position (xyz)
            res = split(tline);
            x = str2num(res{2});
            y = str2num(res{3});
            z = str2num(res{4});
            content.recPosRaw = [x,y,z]';
       elseif linenum == 23
            res = split(tline);
            year = str2num(res{2});Month = str2num(res{3});Day = str2num(res{4});
            Hour = str2num(res{5});Minute = str2num(res{6});Second = str2num(res{7});
            content.date = struct('year',year,'Month',Month,'Day',Day, ...
                                  'Hour',Hour,'Minute',Minute,'Second',Second);            
        elseif linenum >= (shiftbase)
            if (numel(tline) < 15) || strcmp(tline(1:3), 'EOF') == 1
                continue;
            end
            if isempty(abs(str2num(tline(1:15)))>0)
%                 disp(linenum);
%                 disp(tline);
                continue;
            end
%             disp(linenum);
            if ~isempty(find((tline == 'G')==1))
                linenum = shiftbase;
            end
            if linenum == shiftbase
                clear section;
                kk = 1;
                incompleteData = 0;
                Year = (tline(1:3));
                Month = tline(5:6);
                Day = tline(8:9);
                Hour = tline(11:12);
                Minute = tline(14:15);
                Second = tline(17:26);
                satNum = str2num(tline(31:32));
                SNRs = cell(satNum,1);
                for j = 1:satNum
                    id = 33+j*3-3;
                    type = tline(id);
                    num = str2num(tline(id+1:id+2));
                    SNRs{j}=strcat(type,num2str(num));
                end
                section.year = str2num(Year);
                section.Month = str2num(Month);
                section.Day = str2num(Day);
                section.Hour = str2num(Hour);
                section.Minute = str2num(Minute);
                section.Second = str2num(Second);
                section.SNRs = SNRs;
            elseif (mod(linenum-shiftbase,2)==1)
                res = split(tline);
                if size(res,1) ~= 6
                    incompleteData = 1;
                else
                    satPos.L1 = str2num(res{2});
                    satPos.L2 = str2num(res{3});
                    satPos.C1 = str2num(res{4});
                    satPos.P2 = str2num(res{5});
                    satPos.P1 = str2num(res{6});
                end
            else
                res = split(tline);
                if size(res,1) ~= 3
                    incompleteData = 1;
                else
                    satPos.D1 = str2num(res{2});
                    satPos.D2 = str2num(res{3});
                    section.satPos(kk,1) = satPos;
                    kk = kk + 1;
                end
            end
            
            %% valid
            if linenum == (shiftbase+satNum*2)
                shiftbase = shiftbase + satNum*2 + 1;
                if incompleteData == 0
                    content.sections{k,1} = section;
                    k = k + 1;
                end
            end
        end        
        linenum = linenum + 1;
    end
    varargout{1} = content;
end