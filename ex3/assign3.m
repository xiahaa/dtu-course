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

% Last Modified by GUIDE v2.5 22-Sep-2018 17:06:32

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
cla(handles.axes1);
axis(handles.axes1, 'equal');
grid(handles.axes1, 'on');
box(handles.axes1, 'on');  
hold(handles.axes1, 'on')
[x,y,z] = sphere;
x = x.*6371000;
y = y.*6371000;
z = z.*6371000;
obj.HEarth =  gobjects(1);
obj.HEarth = surf(x,y,z,'Parent',handles.axes1); % earth
xlabel(handles.axes1, 'x (m)');
ylabel(handles.axes1, 'y (m)');
zlabel(handles.axes1, 'z (m)');

obj.point =  gobjects(1);
obj.point = plot3(0,0,0,'m*','MarkerSize',10,'Parent',handles.axes1);

obj.HSats = gobjects(100);
obj.HTrajectories = gobjects(100);
obj.texts = gobjects(100);
for i = 1:100
    obj.HSats(i) = plot3(0,0,0,'rd','MarkerSize',10,'Parent',handles.axes1);
    line = [[0, 0, 0];[0,0,0]];
    obj.HTrajectories(i) = plot3(line(:,1),line(:,2),line(:,3),'LineWidth',3,'Parent',handles.axes1);
    obj.texts(i) = text(0,0,0,strcat('dist=',num2str(0)),'Interpreter','latex','Parent',handles.axes1);
end
view(3);
hold(handles.axes1, 'off');    

set(handles.axes1,'UserData', obj);
addpath('../utils/');

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


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    handles.path = uigetfile('*.sp3','Open a SP3 file...');
    disp(handles.path);
    content = sp3fileParser(handles.path);
    display = sprintf('Coordinate System: %s\nAgency: %s\n%s-%s-%s:%s:%s:%s\nNum of Satellites: %s\n', ...
        content.firstline.CoordianteSys,content.firstline.Agency, ...
        content.firstline.YearStart, ...
        content.firstline.MonthStart,content.firstline.DayStart, content.firstline.HoureStart, ...
        content.firstline.MinuteStart,content.firstline.SecondStart, num2str(content.satNum));
    i = 1;
    while i <= content.satNum
        display = sprintf('%sName: %s, Accuracy: %s    Name: %s, Accuracy: %s\n',display, ...
            content.satNames(i,:),content.accuracy(i,:), content.satNames(i+1,:),content.accuracy(i+1,:));
        i = i + 2;
    end
    
    set(handles.text5,'String',display);
    set(handles.text5,'HorizontalAlignment','left');
    
    lists = get(handles.listbox1,'String');
    for i = 1:size(content.sections,1)
        list=sprintf('%d-%d-%d:%d:%d:%f',content.sections{i}.year, ...
            content.sections{i}.Month, content.sections{i}.Day, ...
            content.sections{i}.Hour, content.sections{i}.Minute, ...
            content.sections{i}.Second);
        lists = [lists;cellstr(list)];
    end
    set(handles.listbox1,'String',lists);
    hObject.UserData = content;

% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
    

% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    sp3 = handles.pushbutton1.UserData;
    handles.epoch = get(handles.listbox1,'Value');
    lati = str2double(get(handles.edit1,'String'));
    longi = str2double(get(handles.edit2,'String'));
    alti = str2double(get(handles.edit3,'String'));
    data = zeros(sp3.satNum,3);
    section = sp3.sections{handles.epoch};
    for i = 1:sp3.satNum
        data(i,:) = [section.satPos(i).x,section.satPos(i).y,section.satPos(i).z];
    end
    data = data.*1000;
    consParams = struct('a',6378137.0,'f',1/298.257223563);
    [xo,yo,zo] = llhtoCartesian(lati, longi, alti, consParams);
    
    %% compute
    azimuths = zeros(size(data,1),1);
    zeniths = zeros(size(data,1),1);
    visibles = zeros(size(data,1),1);
    dists = ones(size(data,1),1).*-1;
    thresh = get(handles.slider4,'Value');
    for i = 1:size(data,1)
        satPos = data(i,:);
        dx = satPos(1) - xo;
        dy = satPos(2) - yo;
        dz = satPos(3) - zo;
        [e,n,u] = WGS842ENU(lati, longi, dx, dy, dz);
        
        %% compute azimuth and zenith
        azimuth = rad2deg(atan2(e,n));
        zenith = 90 - rad2deg(asin(u / sqrt(e^2+n^2+u^2)));
        elevation = rad2deg(asin(u / sqrt(e^2+n^2+u^2)));
        azimuths(i) = azimuth;
        zeniths(i) = zenith;
        
        if elevation >= thresh
            % identify as visible
            visibles(i) = 1;
            dists(i) = sqrt(dx^2+dy^2+dz^2);
        else
            visibles(i) = 0;
        end
    end
    obj = handles.axes1.UserData;
    obj.point.XData = xo;obj.point.YData = yo;obj.point.ZData = zo;
    for i = 1:size(data,1)
        if visibles(i) == 1
            obj.HSats(i).XData = data(i,1);
            obj.HSats(i).YData = data(i,2);
            obj.HSats(i).ZData = data(i,3);
            obj.HSats(i).Marker = 'd';
            obj.HSats(i).Color = 'r';

            line = [[xo, yo, zo];[data(i,1),data(i,2),data(i,3)]];
            obj.HTrajectories(i).XData = line(:,1);
            obj.HTrajectories(i).YData = line(:,2);
            obj.HTrajectories(i).ZData = line(:,3);
            
            obj.texts(i).Position = [0.5*(xo+data(i,1)),0.5*(yo+data(i,2)),0.5*(zo+data(i,3))];
            obj.texts(i).String = strcat('dist=',num2str(dists(i)));
        else
            line = [[0, 0, 0];[0,0,0]];
            obj.HTrajectories(i).XData = line(:,1);
            obj.HTrajectories(i).YData = line(:,2);
            obj.HTrajectories(i).ZData = line(:,3);
            obj.texts(i).Position = [0,0,0];
            obj.texts(i).String = strcat('','');
            
            obj.HSats(i).XData = data(i,1);
            obj.HSats(i).YData = data(i,2);
            obj.HSats(i).ZData = data(i,3);
            obj.HSats(i).Marker = 's';
            obj.HSats(i).Color = 'k';
        end
    end

    
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
