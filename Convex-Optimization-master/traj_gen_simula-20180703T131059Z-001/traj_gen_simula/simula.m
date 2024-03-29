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

% Last Modified by GUIDE v2.5 24-Jun-2019 23:00:34

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
global count m keyframe;
count = 0;
m = num;
keyframe = [];

set(gcf, 'WindowButtonMotionFcn', @mouseMove1);
set(gcf, 'WindowButtonDownFcn', @mouseClick1);
i = 0;
while count < num
%     C = get(gca, 'CurrentPoint');
%     title(gca,['(X,Y) = (', num2str(C(1,1),2), ', ',num2str(C(1,2),2), ')']);
%     [x, y] = ginput(1);
    if i~= count
        i = count;
        x = keyframe(1,end);y = keyframe(2,end);
        plot(x,y,'ro','Parent',handles.axes1);grid on;
        text(x,y+1,int2str(i),'FontSize',14, ...
            'Parent',handles.axes1);
        hold(handles.axes1,'on');
    end
    drawnow
%     waypts=[waypts;[x y]];
%     count = count + 1;
end
% draw the end
if i~= count
    i = count;
    x = keyframe(1,end);y = keyframe(2,end);
    plot(x,y,'ro','Parent',handles.axes1);grid on;
    text(x,y+1,int2str(i),'FontSize',14, ...
        'Parent',handles.axes1);
    hold(handles.axes1,'on');
end
drawnow

waypts = keyframe';
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
    
    wayptsz = get(handles.pushbutton4,'UserData');
    initvz = str2double(get(handles.initvz,'String'));
    initaz = str2double(get(handles.initaz,'String'));
    endvz = str2double(get(handles.endvz,'String'));
    endaz = str2double(get(handles.endaz,'String'));
    maxvz = str2double(get(handles.maxvz,'String'));
    maxaz = str2double(get(handles.maxaz,'String'));
    
    l = str2double(get(handles.corridor,'String'));% corridor

    initv = [initvx, initvy];
    endv = [endvx,endvy];
    inita = [initax,initay];
    enda = [endax,enday];

    [polyCoeffs,realt] = cvx_project_traj_gen_solver_refine(waypts,initv,inita,endv,enda,maxv,maxa, ...
        wayptsz,initvz,endvz,initaz,endaz,maxvz,maxaz,l);
    order = 5;
    [pts,vts,ats,tss]=sample_pva(polyCoeffs, realt, order);
    
    trajgenRes.coeff = polyCoeffs;    trajgenRes.time = realt;
    set(hObject,'UserData',trajgenRes);
    
    plot(pts(:,1),pts(:,2),'b-','Parent',handles.axes1);
    set(get(handles.axes1,'title'),'string','2D-Trajecory');

    plot3(pts(:,1),pts(:,2),pts(:,3),'b-','Parent',handles.axes9);hold on;
    plot3(waypts(:,1),waypts(:,2),wayptsz(:,1),'ro','MarkerSize',3,'Parent',handles.axes9);
    set(handles.axes9,'View',[-30.4000 45.2000]);
    set(handles.axes9,'Xlabel',text('String','X'));
    set(handles.axes9,'Ylabel',text('String','Y'));
    set(handles.axes9,'Zlabel',text('String','Z'));
% view([-30.4000 45.2000]);
    set(get(handles.axes9,'title'),'string','3D-Trajecory');
    
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
    %% z
    plot(tss,vts(:,3),'b-','Parent',handles.axes7);
    hold(handles.axes7,'on');
    grid(handles.axes7,'on');
    plot(tss,ones(numel(vts(:,3)),1).*maxvz,'r.-','Parent',handles.axes7);
    plot(tss,ones(numel(vts(:,3)),1).*-maxvz,'r.-','Parent',handles.axes7);
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
    
    %%z
    plot(tss,ats(:,3),'b-','Parent',handles.axes8);
    hold(handles.axes8,'on');
    grid(handles.axes8,'on');
    plot(tss,ones(numel(ats(:,3)),1).*maxaz,'r.-','Parent',handles.axes8);
    plot(tss,ones(numel(ats(:,3)),1).*-maxaz,'r.-','Parent',handles.axes8);
    legend('Generated','Maximum');
    
    figure(1)
    cdata = sqrt(vts(:,1).^2+vts(:,2).^2+vts(:,3).^2);
    %% plot 3D trajectory and generate necessary 3D corridor
    h1 = cline(pts(:,1),pts(:,2), pts(:,3),cdata);hold on; grid on;
    set(h1, 'LineWidth', 3);
%     h2 = plot_3D_corridor_and_compute_distance_to_hyperplane(waypts, wayptsz, l, polyCoeffs);
    [h2, ts, ds1, ds2] = plot_3D_corridor_and_compute_distance_to_hyperplane(waypts, wayptsz, l, polyCoeffs, realt);
%     h1 = plot(pts(:,1),pts(:,2),'b-');hold on; grid on;
%     for i = 1:(num-1)
%         x1 = waypts(i,:);
%         x2 = waypts(i+1,:);
%         [a,b,c] = compute2Dline(x1, x2);
%         [x1l,y1l,x1r,y1r] = computeCorridor(a,b,c,x1(1),x1(2),l);
%         [x2l,y2l,x2r,y2r] = computeCorridor(a,b,c,x2(1),x2(2),l);
% 
%         %% line 1
% %         hcorridor = plot([x1(1) x2(1)], [x1(2) x2(2)], 'blue','LineStyle','-');
% 
%         %% line 2
%         hcorridor = plot([x1l x2l], [y1l y2l], 'red','LineStyle','--');
% 
%         %% line 3
%         plot([x1r x2r], [y1r y2r], 'red','LineStyle','--');
%     end
    title('3D Trajecory','Interpreter','latex');
    xlabel('x: (m)','Interpreter','latex');
    ylabel('y: (m)','Interpreter','latex');
    zlabel('z: (m)','Interpreter','latex');
% %     title('Trajectory','Interpreter','latex');
    legend([h1 h2],'Trajectory', 'Corridor', 'Interpreter', 'latex');
    cv=colorbar;
    cv.Label.String = '$velocity\ (m/s)$';
    cv.Label.Interpreter = 'latex';
    cv.FontSize = 12;
    
    figure(4)
    ds = [ds1 ds2];
    cd = max(ds,[],2);
    mapd = colormap(jet);
%     hd1 = ccplot(ts,ds1,cd,mapd);grid on;hold on;
%     hd2 = ccplot(ts,ds2,cd,mapd);grid on;hold on;
    hd1 = plot(ts,ds1, 'm-','LineWidth',3);grid on;hold on;
    hd2 = plot(ts,ds2, 'b-','LineWidth',3);grid on;hold on;
    hdc = plot(ts,ones(numel(ds2),1).*l,'r--','LineWidth', 3);
    legend([hd1 hd2 hdc], 'Distance to Plane 1', 'Distance to Plane 2', 'Maximum Tolerant Distance');
    xlabel('Time: (s)','Interpreter','latex');
    ylabel('$distance:\ (m)$','Interpreter','latex');
    title('Corridor constraint','Interpreter','latex');
    
    figure
    
    h = cline(pts(:,1),pts(:,2),pts(:,3), cdata);
    
    figure(2)
    vx1 = subplot(2,1,1);
%     c = sqrt(vts(:,1).^2);
%     map = colormap(jet);
    h11 = plot(tss,vts(:,1),'g-','LineWidth',3);grid on;hold on;
    h12 = plot(tss,vts(:,2),'b-','LineWidth',3);grid on;hold on;
    h13 = plot(tss,ones(numel(vts(:,1)),1).*maxv,'r--','LineWidth', 3);
    plot(tss,ones(numel(vts(:,1)),1).*-maxv,'r--','LineWidth', 3);
    legend([h11 h12 h13],'v_x','v_y','Maximum','Interpreter','latex');
    xlabel('Time: (s)','Interpreter','latex');
    ylabel('$velocity: (m/s)$','Interpreter','latex');
    title('Velocity Profile','Interpreter','latex');
    subplot(2,1,2);
%     c = sqrt(vts(:,2).^2);
%     map = colormap(jet);
    h13 = plot(tss,vts(:,3),'c-','LineWidth',3);grid on;hold on;
    h14 = plot(tss,ones(numel(vts(:,2)),1).*maxvz,'r--','LineWidth', 3);
    plot(tss,ones(numel(vts(:,2)),1).*-maxvz,'r--','LineWidth', 3);
    legend([h13 h14],'v_z','Maximum','Interpreter','latex');
    xlabel('Time: (s)','Interpreter','latex');
    ylabel('$velocity: (m/s)$','Interpreter','latex');
%     subplot(3,1,3);
%     c = sqrt(vts(:,3).^2);
%     map = colormap(jet);
%     h15 = ccplot(tss,vts(:,3),c,map);grid on;hold on;
%     h16 = plot(tss,ones(numel(vts(:,3)),1).*maxvz,'r--','LineWidth', 3);
%     plot(tss,ones(numel(vts(:,3)),1).*-maxvz,'r--','LineWidth', 3);
%     legend([h15(1) h16],'Generated','Maximum','Interpreter','latex');
%     xlabel('Time: (s)','Interpreter','latex');
%     ylabel('$v_z: (m/s)$','Interpreter','latex');
%     hp1 = get(subplot(3,1,1),'Position');
%     hp2 = get(subplot(3,1,3),'Position');
%     cv = colorbar('Position', [hp1(1)+hp1(3)+0.01 hp2(2) 0.01 hp2(4)*3.8]);
%     cv.Label.String = '$velocity\ (m/s)$';
%     cv.Label.Interpreter = 'latex';
%     cv.FontSize = 12;
%     caxis([0,max(sqrt(vts(:,2).^2+vts(:,1).^2+vts(:,3).^2))]);
    
    figure    
    vx2 = subplot(2,1,1);
%     ca = sqrt(ats(:,1).^2);
%     map1 = colormap(jet);
    h25 = plot(tss,ats(:,1),'g-', 'LineWidth', 3);grid on;hold on;
    h26 = plot(tss,ats(:,2),'b-', 'LineWidth', 3);grid on;hold on;
    plot(tss,ones(numel(ats(:,1)),1).*maxa,'r--', 'LineWidth', 3);
    h27 = plot(tss,ones(numel(ats(:,1)),1).*-maxa,'r--','LineWidth', 3);
    legend([h25 h26 h27],'a_x','a_y','Maximum','Interpreter','latex');
    xlabel('Time: (s)','Interpreter','latex');
    ylabel('$acceleration: (m/s^2)$','Interpreter','latex');
    title('Acceleration Profile','Interpreter','latex');
    subplot(2,1,2);
%     ca = sqrt(ats(:,2).^2);
    h27 = plot(tss,ats(:,3),'c-', 'LineWidth', 3);grid on;hold on;
    h28 = plot(tss,ones(numel(ats(:,2)),1).*maxaz,'r--','LineWidth', 3);
    plot(tss,ones(numel(ats(:,2)),1).*-maxaz,'r--','LineWidth', 3);
    legend([h27 h28],'a_z','Maximum','Interpreter','latex');
    xlabel('Time: (s)','Interpreter','latex');
    ylabel('$acceleration: (m/s^2)$','Interpreter','latex');
%     subplot(3,1,3);
%     ca = sqrt(ats(:,3).^2);
%     h29 = ccplot(tss,ats(:,3),ca,map1);grid on;hold on;
%     h210 = plot(tss,ones(numel(ats(:,3)),1).*maxa,'r--','LineWidth', 3);
%     plot(tss,ones(numel(ats(:,2)),1).*-maxa,'r--','LineWidth', 3);
%     legend([h29(1) h210],'Generated','Maximum','Interpreter','latex');
%     xlabel('Time: (s)','Interpreter','latex');
%     ylabel('$a_z: (m/s^2)$','Interpreter','latex');
%     
%     hp3 = get(subplot(3,1,1),'Position');
%     hp4 = get(subplot(3,1,3),'Position');
%     cv1 = colorbar('Position', [hp3(1)+hp3(3)+0.01 hp4(2) 0.01 hp4(4)*3.8]);
%     cv1.Label.String = '$acceleration\ (m/s^2)$';
%     cv1.Label.Interpreter = 'latex';
%     cv1.FontSize = 12;
%     caxis([0,max(sqrt(ats(:,1).^2 + ats(:,2).^2 + ats(:,3).^2))]);
    
    

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
cla(handles.axes6);
cla(handles.axes7);
cla(handles.axes8);
cla(handles.axes9);
keyframe = [];
count = 0;

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



function maxvz_Callback(hObject, eventdata, handles)
% hObject    handle to maxvz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxvz as text
%        str2double(get(hObject,'String')) returns contents of maxvz as a double


% --- Executes during object creation, after setting all properties.
function maxvz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxvz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxaz_Callback(hObject, eventdata, handles)
% hObject    handle to maxaz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxaz as text
%        str2double(get(hObject,'String')) returns contents of maxaz as a double


% --- Executes during object creation, after setting all properties.
function maxaz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxaz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function initvz_Callback(hObject, eventdata, handles)
% hObject    handle to initvz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of initvz as text
%        str2double(get(hObject,'String')) returns contents of initvz as a double


% --- Executes during object creation, after setting all properties.
function initvz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initvz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function endvz_Callback(hObject, eventdata, handles)
% hObject    handle to endvz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endvz as text
%        str2double(get(hObject,'String')) returns contents of endvz as a double


% --- Executes during object creation, after setting all properties.
function endvz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to endvz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function endaz_Callback(hObject, eventdata, handles)
% hObject    handle to endaz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endaz as text
%        str2double(get(hObject,'String')) returns contents of endaz as a double


% --- Executes during object creation, after setting all properties.
function endaz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to endaz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function initaz_Callback(hObject, eventdata, handles)
% hObject    handle to initaz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of initaz as text
%        str2double(get(hObject,'String')) returns contents of initaz as a double


% --- Executes during object creation, after setting all properties.
function initaz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initaz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
num = str2double(get(handles.pts,'String'));
wayptsz = [];
hold(handles.axes6,'on');
global count m keyframez;
count = 0;
m = num;
keyframez = [];

set(gcf, 'WindowButtonMotionFcn', @mouseMove2);
set(gcf, 'WindowButtonDownFcn', @mouseClick2);
i = 0;

while count < num
    if i~= count
        i = count;
        z = keyframez(2,end);
        x = keyframez(1,end);
        plot([x-0.1 x+0.1],[z z],'r-','LineWidth',2,'Parent',handles.axes6);
        grid on;
        text(x, z+1, num2str(i), 'FontSize',14, ...
            'Parent',handles.axes6);
        hold(handles.axes6,'on');
    end
    drawnow
end
% draw the end
if i~= count
    i = count;
    z = keyframez(2,end);
    x = keyframez(1,end);
    plot([x-0.1 x+0.1],[z z],'r-','LineWidth',2,'Parent',handles.axes6);
    grid on;
    text(x, z+1, num2str(i), 'FontSize',14, ...
        'Parent',handles.axes6);
    hold(handles.axes6,'on');
end
drawnow
wayptsz = keyframez(2,:)';
set(hObject,'UserData',wayptsz);



function corridor_Callback(hObject, eventdata, handles)
% hObject    handle to corridor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of corridor as text
%        str2double(get(hObject,'String')) returns contents of corridor as a double


% --- Executes during object creation, after setting all properties.
function corridor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to corridor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
waypts = get(handles.pushbutton1,'UserData');
num = size(waypts,1);
l = str2double(get(handles.corridor,'String'));
hold(handles.axes1,'on');
grid(handles.axes1,'on');
for i = 1:(num-1)
    x1 = waypts(i,:);
    x2 = waypts(i+1,:);
    [a,b,c] = compute2Dline(x1, x2);
    [x1l,y1l,x1r,y1r] = computeCorridor(a,b,c,x1(1),x1(2),l);
    [x2l,y2l,x2r,y2r] = computeCorridor(a,b,c,x2(1),x2(2),l);
    
    %% line 1
    plot([x1(1) x2(1)], [x1(2) x2(2)], 'blue','LineStyle','-','Parent',handles.axes1);
    
    %% line 2
    plot([x1l x2l], [y1l y2l], 'red','LineStyle','--','Parent',handles.axes1);
    
    %% line 3
    plot([x1r x2r], [y1r y2r], 'red','LineStyle','--','Parent',handles.axes1);
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    addpath('ctrl') % todo
    h1 = figure;
    res = get(handles.pushbutton2,'UserData');
    
    order = 5;
    [pts,vts,ats,tss]=sample_pva(res.coeff, res.time, order);
    
    plot3(pts(:,1),pts(:,2),pts(:,3),'m-', 'LineWidth', 3);hold on;
    set(gca,'Xlabel',text('String','X'));
    set(gca,'Ylabel',text('String','Y'));
    set(gca,'Zlabel',text('String','Z'));
    set(get(gca,'title'),'string','3D-Trajecory');
    ctrl_entry(h1, res.coeff, res.time);
