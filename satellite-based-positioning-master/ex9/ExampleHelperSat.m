classdef ExampleHelperSat < handle
%Create a simulated sat plot.
% Author xiahaa@space.dtu.dk. 
    properties(SetAccess = public)
        %Trace Pose history of the sat
        Traces = [];
        len = 0;
        sats = [];
        init = 0;
    end
    
    properties(Access = private)  
        %HTrajectory Graphics handle for the trajectory of the sat
        HTrajectories
        
        %HEarth Graphics handle for earth
        HEarth
        
        %HParticles Graphics handle for sats
        HSats
        
        %FigureHandle the handle of the figure 
        FigureHandle
    end
    
    
    methods
        function obj = ExampleHelperSat(numSat, sats)
            %ExampleHelperSat Constructor
            obj.Traces = zeros(numSat,1,3);
            for i=1:numSat
                obj.Traces(i,1,:) = sats(i,:);
            end
            
            obj.len = 1;
            obj.sats = zeros(numSat,3);
           
            obj.FigureHandle = figure('Name', 'Satellite Orbit Gif ');
            % clear the figure
            ax = axes(obj.FigureHandle);
            cla(ax)
%             set(gcf,'Menubar','none','Name','Spheres', ...
%                  'NumberTitle','off','Position',[10,350,400,300], ...
%                  'Color',[0 0 0]);

            % customize the figure
            obj.FigureHandle.Position = [100 100 500 500];
            axis(ax, 'equal');
            grid(ax, 'on');
            box(ax, 'on');         
            
            hold(ax, 'on')
            [x,y,z] = sphere;
            x = x.*6371000;
            y = y.*6371000;
            z = z.*6371000;
           
            obj.HEarth =  surf(x,y,z); % earth
            color = jet(numSat);
            markers = ['o','*','s','d'];
            obj.HSats = gobjects(numSat);
            obj.HTrajectories = gobjects(numSat);
            for i = 1:numSat
                if mod(i,4) == 0
                    id = 4;
                else
                    id = mod(i,4);
                end
                obj.HSats(i) = scatter3(ax, 0,0,0,'MarkerEdgeColor',color(i,:), 'Marker', markers(id));
                obj.HTrajectories(i) = plot3(ax, obj.Traces(i,:,1), obj.Traces(i,:,2), obj.Traces(i,:,3),'Color',color(i,:),'LineWidth',2);
            end
            
            title(ax, strcat('Satellite Orbit Gif', 't = 0'));
            xlabel(ax, 'x (m)');
            ylabel(ax, 'y (m)');
            zlabel(ax, 'z (m)');
            hold(ax, 'off');        
            
            view(3)
        end
        
        function updatePlot(obj, sats, t)
            % updatePlot            
            % render sats
            obj.sats = sats;
            for i = 1:size(sats,1)
                obj.HSats(i).XData = obj.sats(i,1);
                obj.HSats(i).YData = obj.sats(i,2);
                obj.HSats(i).ZData = obj.sats(i,3);
            end
            
            if obj.len == 100
                obj.Traces(:,1:99,:) = obj.Traces(:,2:100,:);
                id = 100;
            else
                id = obj.len + 1;
                obj.len = obj.len + 1;
            end
            for i = 1:size(sats,1)
                obj.Traces(i,id,:) = obj.sats(i,:);
            end
            
            % draw trajectories
            for i = 1:size(sats,1)
               obj.HTrajectories(i).XData =  obj.Traces(i,:,1);
               obj.HTrajectories(i).YData =  obj.Traces(i,:,2);
               obj.HTrajectories(i).ZData =  obj.Traces(i,:,3);
            end
            
            ax = get(obj.FigureHandle, 'currentaxes');
            title(ax, strcat('Satellite Orbit Gif ',' t = ',num2str(t)));
            
%             oh=obj.HEarth;
            %Spins about z axis.
%             axis on;
%             rotate(obj.HEarth,[0 0 1],t/3600*15);
        end
        
    end
    
end