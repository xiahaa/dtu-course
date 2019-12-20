function varargout = assign3(varargin)
% ASSIGN3 MATLAB code for assign3.fig
%      ASSIGN3, by itself, creates a new ASSIGN3 or raises the existing
%      singleton*.
%
%      H = ASSIGN3 returns the handle to a new ASSIGN3 or the handle to
%      the existing singleton*.
%
%      ASSIGN3('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ASSIGN3.M with the given input arguments.
%
%      ASSIGN3('Property','Value',...) creates a new ASSIGN3 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before assign3_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to assign3_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help assign3

% Last Modified by GUIDE v2.5 04-Dec-2018 23:23:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @assign3_OpeningFcn, ...
                   'gui_OutputFcn',  @assign3_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before assign3 is made visible.
function assign3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to assign3 (see VARARGIN)

% Choose default command line output for assign3
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

value = get(handles.slider1,'Value');
set(handles.edit1,'String',num2str(value));
value = get(handles.slider2,'Value');
set(handles.edit2,'String',num2str(value));
value = get(handles.slider3,'Value');
set(handles.edit3,'String',num2str(value));
value = get(handles.slider4,'Value');
set(handles.edit4,'String',num2str(value));
%% h,m,s
value = get(handles.slider5,'Value');
set(handles.edit5,'String',num2str(value));
value = get(handles.slider6,'Value');
set(handles.edit6,'String',num2str(value));
value = get(handles.slider7,'Value');
set(handles.edit7,'String',num2str(value));
%% y,m,d
value = get(handles.slider8,'Value');
set(handles.edit18,'String',num2str(value));
value = get(handles.slider9,'Value');
set(handles.edit19,'String',num2str(value));
value = get(handles.slider10,'Value');
set(handles.edit20,'String',num2str(value));
%% llh
set(handles.edit1,'String',num2str(55.78575300466123));
set(handles.edit2,'String',num2str(12.525384183973078));
set(handles.edit3,'String',num2str(40));
set(handles.slider1,'Value',(55.78575300466123));
set(handles.slider2,'Value',(12.525384183973078));
set(handles.slider3,'Value',(40));
%% ymd
set(handles.edit18,'String',num2str(2018));
set(handles.edit19,'String',num2str(12));
set(handles.edit20,'String',num2str(4));
set(handles.slider8,'Value',(2018));
set(handles.slider9,'Value',(12));
set(handles.slider10,'Value',(4));
%% hms
set(handles.edit5,'String',num2str(12));
set(handles.edit6,'String',num2str(0));
set(handles.edit7,'String',num2str(0));
set(handles.slider5,'Value',(12));
set(handles.slider6,'Value',(0));
set(handles.slider7,'Value',(0));

set(handles.slider11,'Value',(0));


cla(handles.axes1);
% axis(handles.axes1, 'equal');
% grid(handles.axes1, 'on');
% box(handles.axes1, 'on');  
% hold(handles.axes1, 'on')
% [x,y,z] = sphere;
% x = x.*6371000;
% y = y.*6371000;
% z = z.*6371000;
% obj.HEarth =  gobjects(1);
% obj.HEarth = surf(x,y,z,'Parent',handles.axes1); % earth
% xlabel(handles.axes1, 'x (m)');
% ylabel(handles.axes1, 'y (m)');
% zlabel(handles.axes1, 'z (m)');
% obj.point =  gobjects(1);
% obj.point = plot3(0,0,0,'m*','MarkerSize',10,'Parent',handles.axes1);
% obj.HSats = gobjects(100);
% obj.HTrajectories = gobjects(100);
% obj.texts = gobjects(100);
% for i = 1:100
%     obj.HSats(i) = plot3(0,0,0,'rd','MarkerSize',10,'Parent',handles.axes1);
%     line = [[0, 0, 0];[0,0,0]];
%     obj.HTrajectories(i) = plot3(line(:,1),line(:,2),line(:,3),'LineWidth',3,'Parent',handles.axes1);
%     obj.texts(i) = text(0,0,0,strcat('dist=',num2str(0)),'Interpreter','latex','Parent',handles.axes1);
% end
% view(3);
% hold(handles.axes1, 'off');    
% set(handles.axes1,'UserData', obj);
addpath('../utils/');
addpath('../utils/3rdparty');

% UIWAIT makes assign3 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = assign3_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    year = str2double(get(handles.edit18,'String'));
    month = str2double(get(handles.edit19,'String'));
    day = str2double(get(handles.edit20,'String'));
    hour = str2double(get(handles.edit5,'String'));
    minute = str2double(get(handles.edit6,'String'));
    second = str2double(get(handles.edit7,'String'));
    
    %% simulated midnight utc time
    utcMidNight.year = round(year);
    utcMidNight.month = round(month);
    utcMidNight.day = round(day);
    utcMidNight.hour = round(hour); 
    utcMidNight.minute = round(minute);   
    utcMidNight.second = round(second);
    jdMidNight = juliandate(utcMidNight);
    gast = jd2gast(jdMidNight);%% degree
    
    %% position
    lat = str2double(get(handles.edit1,'String'));
    lon = str2double(get(handles.edit2,'String'));
    alti = str2double(get(handles.edit3,'String'));

    consParams = struct('a',6378137.0,'f',1/298.257223563);
    [xo,yo,zo] = llhtoCartesian(lat, lon, alti, consParams);
    rec_xyz = [xo,yo,zo];
    t = get(handles.slider11,'Value');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    theta = gast + t / 3600 * 15;
    R = consR(theta,3);
        
    GalileoKepParams;
    numGalileoSat = size(GalileoKeplerTb,1);    
    %% GPS
    satGPS = {'A1','A2','A3','A4', ...
              'B1','B2','B3','B4', ...
              'C1','C2','C3','C4', ...
              'D1','D2','D3','D4', ...
              'E1','E2','E3','E4', ...
              'F1','F2','F3','F4'}';
    numGPSSat = size(satGPS,1);    
    
    GalileosatPos = zeros(numGalileoSat,3);
    GPSsatPos = zeros(numGPSSat,3);
    
    
    
    %% Galileo
    for j = 1:numGalileoSat
        [~, ~, ~, i1, i2, i3] = calc_sat_pos_with_Kepler(GalileoKeplerTb(j,:), 0, t);
        GalileosatPos(j,1) = i1;
        GalileosatPos(j,2) = i2;
        GalileosatPos(j,3) = i3;
    end
    for j = 1:numGPSSat
        slotID = satGPS{j};
        [~, ~, ~, i1, i2, i3] = calcSatPosition(slotID, 0, t);
        GPSsatPos(j,1) = i1;
        GPSsatPos(j,2) = i2;
        GPSsatPos(j,3) = i3;
    end
    %% transform
    GalileosatPosWGS = (R * GalileosatPos')';
    GPSsatPosWGS = (R * GPSsatPos')';
    %% azi, elv
    dGalileo = GalileosatPosWGS - repmat(rec_xyz,numGalileoSat,1);
    dGPS = GPSsatPosWGS - repmat(rec_xyz,numGPSSat,1);
    [visibilityGalileo, aziGalileo, zenGalileo, elvGalileo] = checkVisibility(numGalileoSat, dGalileo, lat, lon);
    [visibilityGPS, aziGPS, zenGPS, elvGPS] = checkVisibility(numGPSSat, dGPS, lat, lon);
    %% normalization
    normdGalileo = sqrt(dGalileo(:,1).^2+dGalileo(:,2).^2+dGalileo(:,3).^2);
    normdGPS = sqrt(dGPS(:,1).^2+dGPS(:,2).^2+dGPS(:,3).^2);

    idGalileo = (visibilityGalileo == 1);
    idGPS = (visibilityGPS == 1);

    ndGalileo = dGalileo(idGalileo,:) ./ normdGalileo(idGalileo);
    ndGPS = dGPS(idGPS,:) ./ normdGPS(idGPS);

    %% H
    HGalileo = [ndGalileo ones(size(ndGalileo,1),1)];
    HGPS = [ndGPS ones(size(ndGPS,1),1)];

    MGalileo = HGalileo' * HGalileo;%% 4x4
    MGPS = HGPS' * HGPS;%% 4x4

    HIntegration = [ndGalileo ones(size(ndGalileo,1),1) zeros(size(ndGalileo,1),1); ...
                    ndGPS zeros(size(ndGPS,1),1) ones(size(ndGPS,1),1)];%% since we have two independent clock errors.
    MIntegration = HIntegration'*HIntegration;                                  
    %% PDOP
    PDOPGalileo = PDOP_calc(MGalileo);
    PDOPGPS = PDOP_calc(MGPS);
    PDOPIntegration = PDOP_calc(MIntegration);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(handles.edit10,'String',num2str(PDOPGalileo));
    set(handles.edit11,'String',num2str(PDOPIntegration));
    set(handles.edit14,'String',num2str(PDOPGPS));
    
    blue_color = [0 116 186]/255;
    orange_color = [223 80 35]/255;
    
    axes(handles.axes1); 
    azs1 = [aziGalileo(idGalileo)];
    zns1 = [elvGalileo(idGalileo)];
    PRN1 = 1:numGalileoSat;
    PRN1 = PRN1(idGalileo);
    azs2 = [aziGPS(idGPS)];
    zns2 = [elvGPS(idGPS)];
    PRN2 = 1:numGPSSat;
    PRN2 = PRN2(idGPS);
    skyplot_lite(azs1,zns1,PRN1, blue_color, 'o', ...
                 azs2,zns2,PRN2, orange_color, '*');
%     set(hGal,'LineWidth',2);
%     set(hGal,'Color',blue_color);
    

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    value = get(hObject,'Value');
    set(handles.edit1,'String',num2str(value));

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    value = get(hObject,'Value');
    set(handles.edit2,'String',num2str(value));

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    value = get(hObject,'Value');
    set(handles.edit3,'String',num2str(value));

% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
    value = get(hObject,'String');
    value = str2double(value);
    set(handles.slider1,'Value',value);

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
    value = get(hObject,'String');
    value = str2double(value);
    set(handles.slider2,'Value',value);

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
    value = get(hObject,'String');
    value = str2double(value);
    set(handles.slider3,'Value',value);

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    value = get(hObject,'Value');
    set(handles.edit4,'String',num2str(value));
    
% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
    value = get(hObject,'String');
    value = str2double(value);
    set(handles.slider4,'Value',value);

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider5_Callback(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    value = get(hObject,'Value');
    set(handles.edit5,'String',num2str(value));

% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
    value = get(hObject,'String');
    value = str2double(value);
    set(handles.slider5,'Value',value);

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider6_Callback(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    value = get(hObject,'Value');
    set(handles.edit6,'String',num2str(value));

% --- Executes during object creation, after setting all properties.
function slider6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double
    value = get(hObject,'String');
    value = str2double(value);
    set(handles.slider6,'Value',value);

% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider7_Callback(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    value = get(hObject,'Value');
    set(handles.edit7,'String',num2str(value));

% --- Executes during object creation, after setting all properties.
function slider7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
    value = get(hObject,'String');
    value = str2double(value);
    set(handles.slider7,'Value',value);

% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function slider8_Callback(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    value = get(hObject,'Value');
    set(handles.edit18,'String',num2str(value));

% --- Executes during object creation, after setting all properties.
function slider8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double
    value = get(hObject,'String');
    value = str2double(value);
    set(handles.slider8,'Value',value);

% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider9_Callback(hObject, eventdata, handles)
% hObject    handle to slider9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    value = get(hObject,'Value');
    set(handles.edit19,'String',num2str(value));

% --- Executes during object creation, after setting all properties.
function slider9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double
    value = get(hObject,'String');
    value = str2double(value);
    set(handles.slider9,'Value',value);

% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider10_Callback(hObject, eventdata, handles)
% hObject    handle to slider10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    value = get(hObject,'Value');
    set(handles.edit20,'String',num2str(value));

% --- Executes during object creation, after setting all properties.
function slider10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double
    value = get(hObject,'String');
    value = str2double(value);
    set(handles.slider10,'Value',value);

% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    year = str2double(get(handles.edit18,'String'));
    month = str2double(get(handles.edit19,'String'));
    day = str2double(get(handles.edit20,'String'));
    hour = str2double(get(handles.edit5,'String'));
    minute = str2double(get(handles.edit6,'String'));
    second = str2double(get(handles.edit7,'String'));
    
    %% simulated midnight utc time
    utcMidNight.year = round(year);
    utcMidNight.month = round(month);
    utcMidNight.day = round(day);
    utcMidNight.hour = round(hour); 
    utcMidNight.minute = round(minute);   
    utcMidNight.second = round(second);
    jdMidNight = juliandate(utcMidNight);
    gast = jd2gast(jdMidNight);%% degree
    
    %% position
    lat = str2double(get(handles.edit1,'String'));
    lon = str2double(get(handles.edit2,'String'));
    alti = str2double(get(handles.edit3,'String'));

    consParams = struct('a',6378137.0,'f',1/298.257223563);
    [xo,yo,zo] = llhtoCartesian(lat, lon, alti, consParams);
    rec_xyz = [xo,yo,zo];
        
    GalileoKepParams;
    numGalileoSat = size(GalileoKeplerTb,1);    
    %% GPS
    satGPS = {'A1','A2','A3','A4', ...
              'B1','B2','B3','B4', ...
              'C1','C2','C3','C4', ...
              'D1','D2','D3','D4', ...
              'E1','E2','E3','E4', ...
              'F1','F2','F3','F4'}';
    numGPSSat = size(satGPS,1);    
    
    GalileosatPos = zeros(numGalileoSat,3);
    GPSsatPos = zeros(numGPSSat,3);
    
    PDOPs = [];
    times = [];
    for t = 0:60*15:3600*24
        theta = gast + t / 3600 * 15;
        R = consR(theta,3);
        GalileosatPos = zeros(numGalileoSat,3);
        GPSsatPos = zeros(numGPSSat,3);
        %% Galileo
        for j = 1:numGalileoSat
            [~, ~, ~, i1, i2, i3] = calc_sat_pos_with_Kepler(GalileoKeplerTb(j,:), 0, t);
            GalileosatPos(j,1) = i1;
            GalileosatPos(j,2) = i2;
            GalileosatPos(j,3) = i3;
        end
        for j = 1:numGPSSat
            slotID = satGPS{j};
            [~, ~, ~, i1, i2, i3] = calcSatPosition(slotID, 0, t);
            GPSsatPos(j,1) = i1;
            GPSsatPos(j,2) = i2;
            GPSsatPos(j,3) = i3;
        end
        %% transform
        GalileosatPosWGS = (R * GalileosatPos')';
        GPSsatPosWGS = (R * GPSsatPos')';
        %% azi, elv
        dGalileo = GalileosatPosWGS - repmat(rec_xyz,numGalileoSat,1);
        dGPS = GPSsatPosWGS - repmat(rec_xyz,numGPSSat,1);
        visibilityGalileo = checkVisibility(numGalileoSat, dGalileo, lat, lon);
        visibilityGPS = checkVisibility(numGPSSat, dGPS, lat, lon);
        %% normalization
        normdGalileo = sqrt(dGalileo(:,1).^2+dGalileo(:,2).^2+dGalileo(:,3).^2);
        normdGPS = sqrt(dGPS(:,1).^2+dGPS(:,2).^2+dGPS(:,3).^2);
        
        idGalileo = (visibilityGalileo == 1);
        idGPS = (visibilityGPS == 1);
        
        ndGalileo = dGalileo(idGalileo,:) ./ normdGalileo(idGalileo);
        ndGPS = dGPS(idGPS,:) ./ normdGPS(idGPS);
        
        %% H
        HGalileo = [ndGalileo ones(size(ndGalileo,1),1)];
        HGPS = [ndGPS ones(size(ndGPS,1),1)];
        
        MGalileo = HGalileo' * HGalileo;%% 4x4
        MGPS = HGPS' * HGPS;%% 4x4
        
        HIntegration = [ndGalileo ones(size(ndGalileo,1),1) zeros(size(ndGalileo,1),1); ...
                        ndGPS zeros(size(ndGPS,1),1) ones(size(ndGPS,1),1)];%% since we have two independent clock errors.
        MIntegration = HIntegration'*HIntegration;
        %% PDOP
        PDOPGalileo = PDOP_calc(MGalileo);
        PDOPGPS = PDOP_calc(MGPS);
        PDOPIntegration = PDOP_calc(MIntegration);
        
        %% save
        times = [times;t];
        PDOPs = [PDOPs;[PDOPGalileo PDOPGPS PDOPIntegration]];
    end
    red_color = [153 0 0]/255;
    blue_color = [0 116 186]/255;
    orange_color = [223 80 35]/255;
    font_size = 16;
    axes(handles.axes2); 
    plot(times,PDOPs(:,1),'LineStyle','--','LineWidth',2.5, 'Color', blue_color);grid on;hold on;
    plot(times,PDOPs(:,2),'LineStyle','-.','LineWidth',2.5, 'Color', orange_color);
    plot(times,PDOPs(:,3),'LineStyle','-','LineWidth',2.5, 'Color', red_color);
    title('PDOP over 24 hours','Interpreter','latex');
    xlabel('time: (s)','Interpreter','latex');
    ylabel('PDOP','Interpreter','latex');
    lgnd = legend({'Galileo Only','GPS Only','GPS+Galileo'}, 'Location', 'NorthWest');
    set(lgnd, 'Interpreter', 'Latex','FontSize', font_size);
    colormap summer
    set(gca,'TickLabelInterpreter','latex');
    xlabel('time: (s)','FontSize', font_size, 'Interpreter', 'latex');
    ylabel('PDOP','FontSize', font_size, 'Interpreter', 'latex');
    hold off;   
    
% --- Executes on slider movement.
function slider11_Callback(hObject, eventdata, handles)
% hObject    handle to slider11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
