function path = rrt_star(map3D, sp, ep, maxIter, sizeX, sizeY, sizeZ, resolution)
    %% use rrt_star to find path from sp to ep,
    %% if fail, then path will be empty
    %% else, path will include all waypts
    %% range should be in integer, 
    
    %% test
%     clc;clear all;close all;
%     sizeX = 10;
%     sizeY = 10;
%     sizeZ = 1;
%     resolution = 0.1;
%     sp = [0 0 0];
%     ep = [4 4 0.5];
%     maxIter = 5000;
%     map = [];
%     for i = 1:0.1:2
%         for j = 1:0.1:2
%             for k = 0:0.1:1
%                 map = [map;[i j k]];
%             end
%         end
%     end
%     map3D = robotics.OccupancyMap3D(9.9,'FreeThreshold',0.2,'OccupiedThreshold',0.65);
%     pose = [ 0 0 0 1 0 0 0];
%     maxRange = sizeX * 0.8;
%     insertPointCloud(map3D,pose,map,maxRange);
%     setOccupancy(map3D,map,0.9);
%     inflate(map3D,0.5);
%     show(map3D);hold on;
%     xlim([-5 5]); ylim([-5 5]); zlim([0 2]);
%     grid on;
%     plot3(sp(1),sp(2),sp(3),'ro','MarkerSize',5);
%     plot3(ep(1),ep(2),ep(3),'ro','MarkerSize',5);
%     plot3(map(:,1),map(:,2),map(:,3),'*');

    dims = numel(sp);
    path = [];
    
    scale = 1/resolution;
    sizeX = sizeX * scale;
	sizeY = sizeY * scale;
	sizeZ = sizeZ * scale;
    
    quasiRandomMap = 1:(sizeX*sizeY*sizeZ);

    nodes = [sp];
    costs_from_root = [0];
    node = struct('parent',[],'childrens',[],'id',1);
    rrt_tree = [node];
    
    solExt = true;
    %% rrt_star entry
    for i = 1:1:maxIter
        %% first, check if rrt_star should stop
        if isGoalReached(ep, nodes)
            break;
        end
        %% random pick one point,
        [ranpt,quasiRandomMap,succ] = randomSampling(quasiRandomMap, sizeX, sizeY, sizeZ, scale, map3D, ep);
        if succ
            %% find nn with minimum cost
            [parent, idRadius] = findNN(rrt_tree, ranpt, nodes, costs_from_root);
            %% connect
%             [rrt_tree, newnode, cost_from_root_to_newpt, succ] = connect(parent, ranpt, rrt_tree, map);
            [rrt_tree, newpt, newid, cost_from_root_to_newpt, succ] = connect(parent, ranpt, rrt_tree, map3D, nodes, costs_from_root);
            if succ
                nodes = [nodes;newpt];
                costs_from_root = [costs_from_root; cost_from_root_to_newpt];
                %% rewire
                for j = 1:numel(idRadius)
                    nn = nodes(idRadius(j),:);
                    cost_from_newpt_to_nn = norm(newpt-nn);
                    
                    collision = collisionCheck(newpt, nn, map3D);
                    if collision
                        continue;
                    end
                    
                    if (cost_from_root_to_newpt + cost_from_newpt_to_nn) < costs_from_root(idRadius(j))
                        costs_from_root(idRadius(j)) = cost_from_root_to_newpt + cost_from_newpt_to_nn;
                        %% relink to newpt
                        % delete from parent's children list
                        idp = rrt_tree(idRadius(j)).parent;
                        idc = find(rrt_tree(idp).childrens == idRadius(j));
                        rrt_tree(idp).childrens(idc) = [];
                        % update parent
                        rrt_tree(idRadius(j)).parent = newid;
                        %% iterate to update its childrens
                        costs_from_root = rewire(rrt_tree, idRadius(j), nodes, costs_from_root);
                    end
                end
            end
        else
            break;
        end
    end
    
    %% verify is solution exists
    if isGoalReached(ep, nodes)
        solExt = true;
        
        %% trace path
        pts = nodes;
        dist = repmat(ep,size(pts,1),1) - pts;
        dist = dist(:,1).^2+dist(:,2).^2+dist(:,3).^2;
        [minval,minid] = min(dist);
        pathid = [minid];
        while 1
            idp = rrt_tree(pathid(end)).parent;
            pt = nodes(idp,:);
            dist1 = norm(pt-sp);
            pathid = [pathid;idp];
            if dist1 < 1e-3
                %% found
                break;
            end
        end
        pathid = flipud(pathid);
        %% line-of-sight tracing
        pathlst = [pathid(1)];
        i = 2;
        while i <= numel(pathid)
            pt1 = nodes(pathlst(end),:);
            pt2 = nodes(pathid(i),:);
            collision = collisionCheck(pt1, pt2, map3D);
            if collision
                if isempty(find(pathlst == pathid(i-1)))
                    pathlst = [pathlst;pathid(i-1)];
                else
                    i = i + 1;
                end
            else
                 i = i + 1;
            end
        end
        pathlst = [pathlst;pathid(end)];
        %% assign outputs
        for i = 1:numel(pathlst)
            pt = nodes(pathlst(i),:);
            path = [path;pt];
        end
    else
        solExt = false;
        path = [];
    end
    
%     plot3(path(:,1),path(:,2),path(:,3),'g-');
    
end

function collision = collisionCheck(pt1, pt2, map)
    dist = pt2 - pt1;
    distLen = norm(dist);
    dir = dist / distLen;
    step = 0.1;
    maxIter = floor(distLen / step);
    collisionRadius = 0.75;
    i = 1;
    while i <= maxIter
        ptc = pt1 + dir*i*step;
        %% check if collision 
%         dist = repmat(ptc,size(map,1),1) - map;
%         dist = dist(:,1).^2+dist(:,2).^2+dist(:,3).^2;
%         isCollision = ~isempty(find(dist < collisionRadius^2));
        isCollision = checkOccupancy(map,ptc);
        if isCollision ~= 1
            i = i + 1;
        else
            break;
        end
    end
    if i > maxIter
        collision = false;
    else
        collision = true;
    end
end

function [costs_from_root] = rewire(rrt_tree, id, nodes, costs_from_root)
    if isempty(rrt_tree(id).childrens)
        return;
    else
        new_root_cost = costs_from_root(id);
        for i=1:numel(rrt_tree(id).childrens)
            cost_from_parent_to_children = norm(nodes(rrt_tree(id).childrens(i),:) - nodes(id,:));
            costs_from_root(rrt_tree(id).childrens(i)) = new_root_cost + cost_from_parent_to_children;
            costs_from_root = rewire(rrt_tree, rrt_tree(id).childrens(i), nodes, costs_from_root);
        end
    end
end

function [rrt_tree, newpt, newid, cost_from_root_to_newpt, succ] = connect(parent, ranpt, rrt_tree, map, nodes, costs)
    maxStep = 1;
    ptparent = nodes(parent,:);
    dist = ranpt - ptparent;
    distLen = norm(dist);
    collisionRadius = 0.75;
    stepLen = 0.1;
    
    s = 1;
    dir = dist / distLen;
    while (s*stepLen <= maxStep) && (s*stepLen <= distLen)
        pt = ptparent + dir * s * stepLen;
        %% check if collision 
        isCollision = checkOccupancy(map,pt);
%         dist = repmat(pt,size(map,1),1) - map;
%         dist = dist(:,1).^2+dist(:,2).^2+dist(:,3).^2;
%         isCollision = ~isempty(find(dist < collisionRadius^2));
        if isCollision ~= 1
            s = s + 1;
        else
            s = s - 1;
            break;
        end
    end
    if s > 0
        %% new pt
        newpt = ptparent + dir * s * stepLen;
        %% cost from parent to current
        cost_from_parent_to_newpt = norm(newpt-ptparent);
        %% cost from root to current
        cost_from_root_to_newpt = costs(parent) + cost_from_parent_to_newpt;
        succ = true;
        
        %% link to parent
        id = size(rrt_tree,1) + 1;
        rrt_tree(parent).childrens = [ rrt_tree(parent).childrens; id];
        %% new
        newnode = struct('parent',[parent],'childrens',[],'id',id);
        newid = id;
        rrt_tree = [rrt_tree;newnode];
    else
        newpt = [];
        newid = [];
        cost_from_root_to_newpt = 0;
        succ = false;
    end
end

function reached = isGoalReached(ep, nodes)
    pts = nodes;
    if isempty(nodes)
        reached = false;
    else
        dist = repmat(ep,size(pts,1),1) - pts;
        dist = dist(:,1).^2+dist(:,2).^2+dist(:,3).^2;
        [minval,minid] = min(dist);
        if minval < 0.1
            reached = true;
        else
            reached = false;
        end
    end
end

function [ranpt,quasiRandomMap,succ] = randomSampling(quasiRandomMap, sizeX, sizeY, sizeZ, scale, map, ep)
    prob = rand(1);
    tossGoal = 0.2;
    collisionRadius = 0.75;
    if prob <= tossGoal
        ranpt = ep;
        succ = true;
    else
        while 1
            if isempty(quasiRandomMap)
                ranpt = [];
                succ = false;
                break;
            end
            %% random sampling on grid
            id = randi(size(quasiRandomMap,2));
            quasiRandomMap(id) = [];

            ix = floor(id / (sizeY*sizeZ));
            id = id - ix * (sizeY*sizeZ);
            iy = floor(id / (sizeZ));
            iz = id - iy*sizeZ;

            %% grid to real metrics
            pt = [ix/scale-(sizeX)/(2*scale),iy/scale-(sizeY)/(2*scale),iz/scale];
%             pt = [1.2 1.2 0.3];
            isCollision = checkOccupancy(map,pt);

            %% check if collision 
%             dist = repmat(pt,size(map,1),1) - map;
%             dist = dist(:,1).^2+dist(:,2).^2+dist(:,3).^2;
%             isCollision = ~isempty(find(dist < collisionRadius^2));
% 
            if isCollision ~= 1
                ranpt = pt;
                succ = true;
                break;
            end
        end
    end
end

function [parent, idRadius] = findNN(rrt_tree, ranpt, nodes, costs)
    %% constant
    radius = 1;
    %% pts and costs
    pts = nodes;
    costs_from_root = costs;
    %% nn
    dist = repmat(ranpt,size(pts,1),1) - pts;
    dist = dist(:,1).^2+dist(:,2).^2+dist(:,3).^2;
    idRadius = find(dist < radius^2);
    
    if isempty(idRadius)
        %% no solution in radius, then find nearest
        [minval, idnearest] = min(dist);
        parent = rrt_tree(idnearest).id;
        idRadius = [];
    else
        dbstop if error
        costRadius = costs_from_root(idRadius);
        [val, id] = min(costRadius);
%         disp('idRadius : ');idRadius(id)
%         disp('id : ');id
%         disp('size rrt_tree : '); size(rrt_tree,1)
        parent = rrt_tree(idRadius(id)).id;
%         idRadius = rrt_tree(idRadius).id;
    end
end

 