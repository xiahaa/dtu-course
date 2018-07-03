function varargout = simula(varargin)
% SIMULA MATLAB code for simula.fig
%      SIMULA, by itself, creates a new SIMULA or raises the existing
%      singleton*.
%
%      H = SIMULA returns the handle to a new SIMULA or the handle to
%      the existing singleton*.
%
%      SIMULA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SIMULA.M with the given input arguments.
%
%      SIMULA('Property','Value',...) creates a new SIMULA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before simula_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to simula_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help simula

% Last Modified by GUIDE v2.5 03-Jul-2018 16:10:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @simula_OpeningFcn, ...
                   'gui_OutputFcn',  @simula_OutputFcn, ...
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


% --- Executes just before simula is made visible.
function simula_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to simula (see VARARGIN)
% Choose default command line output for simula
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes simula wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = simula_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
xlim(hObject,[-10, 10]);
ylim(hObject,[-10, 10]);
grid(hObject, 'on');
% hold(handles.axes1,'on');
% Hint: place code in OpeningFcn to populate axes1


% --- Executes on button press in pushbutton1.
function handles = pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
num = str2double(get(handles.pts,'String'));
waypts = [];
hold(handles.axes1,'on');
for i = 1:num
    [x, y] = ginput(1);
    plot(x,y,'ro','Parent',handles.axes1);grid on;
    text(x,y+1,int2str(i),'FontSize',14, ...
        'Parent',handles.axes1);
    hold(handles.axes1,'on');
    waypts=[waypts;[x y]];
end
set(hObject,'UserData',waypts);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    waypts = get(handles.pushbutton1,'UserData');
    num = size(waypts,1);
    initvx = str2double(get(handles.initv,'String'));
    initax = str2double(get(handles.inita,'String'));
    endvx = str2double(get(handles.endv,'String'));
    endax = str2double(get(handles.enda,'String'));
    initvy = str2double(get(handles.initvy,'String'));
    initay = str2double(get(handles.initay,'String'));
    endvy = str2double(get(handles.endvy,'String'));
    enday = str2double(get(handles.enday,'String'));
    maxv = str2double(get(handles.maxv,'String'));
    maxa = str2double(get(handles.maxa,'String'));

    initv = [initvx, initvy];
    endv = [endvx,endvy];
    inita = [initax,initay];
    enda = [endax,enday];

    [polyCoeffs,realt] = cvx_project_traj_gen_solver_refine(waypts,initv,inita,endv,enda,maxv,maxa);
    order = 5;
    [pts,vts,ats,tss]=sample_pva(polyCoeffs, realt, order);

    plot(pts(:,1),pts(:,2),'b-','Parent',handles.axes1);
    set(get(handles.axes1,'title'),'string','Trajecory');
%     title('Trajectory');

    pax11 = handles.axes2;
    pax12 = handles.axes3;
    plot(tss,vts(:,1),'b-','Parent',pax11);
    hold(handles.axes2,'on');
    grid(handles.axes2,'on');
    plot(tss,ones(numel(vts(:,1)),1).*maxv,'r.-','Parent',pax11);
    plot(tss,ones(numel(vts(:,1)),1).*-maxv,'r.-','Parent',pax11);
    legend('Generated','Maximum');
%     set(get(handles.axes2,'title'),'string','Velocity');
    plot(tss,vts(:,2),'b-','Parent',pax12);
    hold(handles.axes3,'on');
    grid(handles.axes3,'on');
    plot(tss,ones(numel(vts(:,2)),1).*maxv,'r.-','Parent',pax12);
    plot(tss,ones(numel(vts(:,2)),1).*-maxv,'r.-','Parent',pax12);
    legend('Generated','Maximum');
    
    pax13 = handles.axes4;
    pax14 = handles.axes5;
    plot(tss,ats(:,1),'b-','Parent',pax13);
    hold(handles.axes4,'on');
    grid(handles.axes4,'on');
    plot(tss,ones(numel(ats(:,1)),1).*maxa,'r.-','Parent',pax13);
    plot(tss,ones(numel(ats(:,1)),1).*-maxa,'r.-','Parent',pax13);
    legend('Generated','Maximum');
%     set(get(handles.axes1,'title'),'string','Acceleration');
    plot(tss,ats(:,2),'b-','Parent',pax14);
    hold(handles.axes5,'on');
    grid(handles.axes5,'on');
    plot(tss,ones(numel(ats(:,2)),1).*maxa,'r.-','Parent',pax14);
    plot(tss,ones(numel(ats(:,2)),1).*-maxa,'r.-','Parent',pax14);
    legend('Generated','Maximum');

function pts_Callback(hObject, eventdata, handles)
% hObject    handle to pts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pts as text
%        str2double(get(hObject,'String')) returns contents of pts as a double


% --- Executes during object creation, after setting all properties.
function pts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function initv_Callback(hObject, eventdata, handles)
% hObject    handle to initv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of initv as text
%        str2double(get(hObject,'String')) returns contents of initv as a double


% --- Executes during object creation, after setting all properties.
function initv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function inita_Callback(hObject, eventdata, handles)
% hObject    handle to inita (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inita as text
%        str2double(get(hObject,'String')) returns contents of inita as a double


% --- Executes during object creation, after setting all properties.
function inita_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inita (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function endv_Callback(hObject, eventdata, handles)
% hObject    handle to endv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endv as text
%        str2double(get(hObject,'String')) returns contents of endv as a double


% --- Executes during object creation, after setting all properties.
function endv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to endv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function enda_Callback(hObject, eventdata, handles)
% hObject    handle to enda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of enda as text
%        str2double(get(hObject,'String')) returns contents of enda as a double


% --- Executes during object creation, after setting all properties.
function enda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to enda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxv_Callback(hObject, eventdata, handles)
% hObject    handle to maxv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxv as text
%        str2double(get(hObject,'String')) returns contents of maxv as a double


% --- Executes during object creation, after setting all properties.
function maxv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxa_Callback(hObject, eventdata, handles)
% hObject    handle to maxa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxa as text
%        str2double(get(hObject,'String')) returns contents of maxa as a double


% --- Executes during object creation, after setting all properties.
function maxa_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxa (see GCBO)
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
cla(handles.axes1);
cla(handles.axes2);
cla(handles.axes3);
cla(handles.axes4);
cla(handles.axes5);

function endvy_Callback(hObject, eventdata, handles)
% hObject    handle to endvy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endvy as text
%        str2double(get(hObject,'String')) returns contents of endvy as a double


% --- Executes during object creation, after setting all properties.
function endvy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to endvy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function enday_Callback(hObject, eventdata, handles)
% hObject    handle to enday (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of enday as text
%        str2double(get(hObject,'String')) returns contents of enday as a double


% --- Executes during object creation, after setting all properties.
function enday_CreateFcn(hObject, eventdata, handles)
% hObject    handle to enday (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function initvy_Callback(hObject, eventdata, handles)
% hObject    handle to initvy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of initvy as text
%        str2double(get(hObject,'String')) returns contents of initvy as a double


% --- Executes during object creation, after setting all properties.
function initvy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initvy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function initay_Callback(hObject, eventdata, handles)
% hObject    handle to initay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of initay as text
%        str2double(get(hObject,'String')) returns contents of initay as a double


% --- Executes during object creation, after setting all properties.
function initay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
