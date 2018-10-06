function content = sp3fileParser(path)
%sp3fileParser: parsing sp3 file.
%   sp3fileParser(path) parses a given sp3 file.
%   path: path to the sp3 file.
%   content: struct
%       firstline: struct
%               versionSym: see sp3 definition.
%               PosOrVel
%               YearStart
%               MonthStart
%               DayStart
%               HoureStart
%               MinuteStart
%               SecondStart
%               NumOfEpochs
%               DataUsed
%               CoordianteSys
%               OrbitType
%               Agency
%       satNum: int, number of satellites.
%       satNames: array, satellite name.
%       accuracy: array, satellite accuracy.
%       sections: cells of struct, each element is a struct defined as
%       follows.
%            section: struct.
%               year
%               Month
%               Day
%               Hour
%               Minute
%               Second
%               satPos: array of satNumx4 [x,y,z,clocktime].
%   Author: xiahaa@space.dtu.dk

    %% fopen file
    fid = fopen(path,"r");
    linenum = 1;
    satNum = 0;
    shiftbase = 23;

    sections = {};
    k = 1;
    while ~feof(fid)
        tline = fgetl(fid);
        if linenum == 1
            firstline = parseFirstLine(tline);
            content.firstline = firstline;
        elseif linenum == 3
            thirdline = parseThirdLine(tline);
            satNum = str2num(thirdline.NumSats);
            content.satNum = satNum;
            content.satNames = thirdline.SatNames;
        elseif linenum >= 4 && linenum <= 7
            res = parse4to7Line(tline);
            content.satNames = [content.satNames;res.SatNames];
        elseif linenum >= 8 && linenum <=12
            res = parse8to12FifthLine(tline);
            if linenum == 8
                content.accuracy = res.Accuracy;
            else
                content.accuracy = [content.accuracy;res.Accuracy];
            end
        elseif linenum >= (shiftbase) && linenum <= (shiftbase+satNum)
            if strcmp(tline(1:3), 'EOF') == 1
                break;
            end
            
            if linenum == shiftbase
                Year = tline(4:7);
                Month = tline(9:10);
                Day = tline(12:13);
                Hour = tline(15:16);
                Minute = tline(18:19);
                Second = tline(21:31);
                section.year = str2num(Year);
                section.Month = str2num(Month);
                section.Day = str2num(Day);
                section.Hour = str2num(Hour);
                section.Minute = str2num(Minute);
                section.Second = str2num(Second);
            else
                satPos = parseSatPos(tline);
                satPos.x = str2num(satPos.x);
                satPos.y = str2num(satPos.y);
                satPos.z = str2num(satPos.z);
                satPos.clock = str2num(satPos.clock);
                section.satPos(linenum - (shiftbase),1) = satPos;
            end
            
            %% valid
            if linenum == (shiftbase+satNum)
                shiftbase = shiftbase + satNum + 1;
                content.sections{k,1} = section;
                k = k + 1;
            end
        end        
        linenum = linenum + 1;
    end
    
end

function satPos = parseSatPos(line)
    satPos.symbol = line(1);
    satPos.name = line(2:4);
    satPos.x = line(5:18);%km
    satPos.y = line(19:32);%km
    satPos.z = line(33:46);%km
    satPos.clock = line(47:60);%micro second
end

function firstline = parseFirstLine(line)
    firstline.versionSym = line(1:2);
    firstline.PosOrVel = line(3);
    firstline.YearStart = line(4:7);
    firstline.MonthStart = line(9:10);
    firstline.DayStart = line(12:13);
    firstline.HoureStart = line(15:16);
    firstline.MinuteStart = line(18:19);
    firstline.SecondStart = line(21:31);
    firstline.NumOfEpochs = line(33:39);
    firstline.DataUsed = line(41:45);
    firstline.CoordianteSys = line(47:51);
    firstline.OrbitType = line(53:55);
    firstline.Agency = line(57:60);
end

function thirdline = parseThirdLine(line)
    thirdline.NumSats = line(5:6);
    thirdline.SatNames = [];
    for i = 1:17
        thirdline.SatNames = [thirdline.SatNames;line(10+(i-1)*3:10+(i-1)*3+2)];
    end
end

function res = parse4to7Line(line)
    res.SatNames = [];
    for i = 1:17
        res.SatNames = [res.SatNames;line(10+(i-1)*3:10+(i-1)*3+2)];
    end
end

function res = parse8to12FifthLine(line)
    res.Accuracy = [];
    for i = 1:17
        res.Accuracy = [res.Accuracy;line(10+(i-1)*3:10+(i-1)*3+2)];
    end
end