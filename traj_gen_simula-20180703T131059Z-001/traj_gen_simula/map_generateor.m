function [output, pcl] = map_generateor(sizeX, sizeY, sizeZ, seed, scale)
    clc;clear all;close all;
    
    %% call c++ interface
    seed = 511;
    sizeX = 14;
    sizeY = 14;
    sizeZ = 2;
    resolution = 0.1;%% 1 grid = 0.25m
    scale = 1/resolution;
    type = 2;
    
    %% obstacles, w_l = minimium width of obstacles in real metrics
    ObsNum = 50;
    w_l = 0.6;
    w_h = 1.5;
    
    complexity = 0.05;
    fill = 0.06;
    fractal = 1;
    attenuation = 0.1;
    
    if type == 1 %% perlin3d
        [px,py,pz] = mexMapGen(seed, sizeX, sizeY, sizeZ, resolution, type, complexity, fill, fractal, attenuation);
    elseif type == 2 %% random
        [px,py,pz] = mexMapGen(seed, sizeX, sizeY, sizeZ, resolution, type,w_l, w_h, ObsNum);
    end
    
%     pclXYZ = randomMapGen(sizeX, sizeY, sizeZ, scale, ObsNum, w_l, w_h);
%     pclXYZ = perlin3D(sizeX, sizeY, sizeZ, resolution);
    
    map3D = robotics.OccupancyMap3D(9.9,'FreeThreshold',0.2,'OccupiedThreshold',0.65);
    pose = [ 0 0 0 1 0 0 0];
    maxRange = sizeX * 0.8;
    pclXYZ = [px py pz];
%     pcshow(pclXYZ,'MarkerSize',1000);
    insertPointCloud(map3D,pose,pclXYZ,maxRange);
    setOccupancy(map3D,pclXYZ,0.9);
    
    map3DInf = robotics.OccupancyMap3D(9.9,'FreeThreshold',0.2,'OccupiedThreshold',0.65);
%     pcshow(pclXYZ,'MarkerSize',1000);
    insertPointCloud(map3DInf,pose,pclXYZ,maxRange);
    setOccupancy(map3DInf,pclXYZ,0.9);
    inflate(map3DInf,0.5);
    
    show(map3D);hold on;
    xlim([-sizeX sizeX]); ylim([-sizeY sizeY]); zlim([0 sizeZ]);
    grid on;
    
    sp = [-7 -7 0];
%     while 1
%         sp(1,1:2) = rand(1,2)*4 - 2;
%         sp(1,3) = 0;
%         isCollision = checkOccupancy(map3D,sp);
%         if isCollision ~= 1
%             break;
%         end
%     end
%     ep = zeros(1,3);
%     while 1
%         ep(1,1:2) = rand(1,2)*4 + 5;
%         ep(1,3) = 1;
%         isCollision = checkOccupancy(map3D,ep);
%         if isCollision ~= 1
%             break;
%         end
%     end
    ep = [8 8 1];
    
%     plot3(sp(1),sp(2),sp(3),'ro','MarkerSize',5);
%     plot3(ep(1),ep(2),ep(3),'bo','MarkerSize',5);
    maxIter = 5000;
    path = rrt_star(map3DInf, sp, ep, maxIter, sizeX, sizeY, sizeZ, resolution);
    
    if ~isempty(path)
        plot3(path(:,1),path(:,2),path(:,3),'g-', 'LineWidth', 3);hold on;
        plot3(path(1:end,1),path(1:end,2),path(1:end,3),'bo');
        initv = [0 0];
        inita = [0 0];
        endv = [0 0];
        enda = [0 0];
        maxv = 10;
        maxa = 6;
        p = path(:,1:2);
        pz = path(:,3);
        initvz = 0;
        endvz = 0;
        initaz = 0;
        endaz = 0;
        maxvz = 3;
        maxaz = 3;
        l = 0.4;
        
        [polyCoeffs,realt]=cvx_project_traj_gen_solver_refine(p,initv,inita,endv,enda, ...
            maxv,maxa,pz,initvz,endvz,initaz,endaz,maxvz,maxaz,l);
        
        order = 5;
        [pts,vts,ats,tss]=sample_pva(polyCoeffs, realt, order);
        
        plot3(pts(:,1),pts(:,2),pts(:,3),'r-','LineWidth',5);hold on;
        title('Perlin Map Test');
        legend('Obstacles', 'RRT* Path','RRT* Way Point', 'Trajectory');
        axis equal
    end


end


function pclXYZ = randomMapGen(sizeX, sizeY, sizeZ, scale, ObsNum, w_l, w_h)
    resolution = 1 / scale;
    
    x_l = -sizeX / (2 * scale);
    x_h =  sizeX / (2 * scale);
    
    y_l = -sizeY / (2 * scale);
    y_h =  sizeY / (2 * scale);

    h_l = 0;
    h_h = sizeZ / scale;
    
%     w_l = 0.6;
%     w_h = 1.5;
    
    pclXYZ = [];
    
    for i = 1:1:ObsNum
        x = rand(1)*(x_h-x_l)+x_l;
        y = rand(1)*(y_h-y_l)+y_l;
        
        w = rand(1)*(w_h-w_l)+w_l;
        h = rand(1)*(h_h-h_l)+h_l;
        
        widNum = ceil(w / resolution);
        heiNum = ceil(h / resolution);
        
        rl = -widNum / 2;
        rh = widNum / 2;
        sl = -widNum / 2;
        sh = widNum / 2;
        
        for r = rl:rh
            for s = sl:sh
                for t = 0:heiNum
                    if (r - rl) * (r - rh + 1) * (s - sl) * (s - sh + 1) * t * (t - heiNum + 1) ==0
                        xs = x + r * resolution;
                        ys = y + s * resolution;
                        zs = t * resolution;
                        pclXYZ = [pclXYZ;[xs ys zs]];
                    end
                end
            end
        end
    end
end

function points = perlin3D(sizeX, sizeY, sizeZ, scale)
    scale = 1 / scale;
	sizeX = sizeX * scale;
	sizeY = sizeY * scale;
	sizeZ = sizeZ * scale;

    complexity = 0.05;
    fill = 0.12;
    fractal = 1;
    attenuation = 0.2;
    
    width  = sizeX * sizeY * sizeZ;
    height = 1;
    points = [];
    
    noise = genPerlinNoise('f');
    
    v = [];
    for (i = 1:sizeX)
        for (j = 1:sizeY)
            for (k = 0:sizeZ-1)
                tnoise = 0;
                for (it = 1:fractal)
                    dfv = 2^it;
                    ta  = attenuation / it;
                    tnoise = tnoise + ta * perlin(noise, dfv * i * complexity, dfv * j * complexity, dfv * k * complexity);
                end
                v = [v;tnoise];
            end
        end
    end
    vs = sort(v);
    tpos = width * (1 - fill);
    tmp  = vs(tpos+1);
    for i = 1:sizeX
        for j = 1:sizeY
            for k = 0:sizeZ-1
                tnoise = 0;
                for (it = 1:fractal)
                    dfv = 2^it;
                    ta  = attenuation / it;
                    tnoise = tnoise + ta * perlin(noise, dfv * i * complexity, dfv * j * complexity, dfv * k * complexity);
                end
                if (tnoise > tmp)
                    points = [points;[(i / scale - sizeX / (2 * scale)) (j / scale - sizeY / (2 * scale)) k/scale]];
                end
            end
        end
    end
end

function PerlinNoise = genPerlinNoise(type)
    if type == 'f'
        PerlinNoise = [151, 160, 137, 91,  90,  15,  131, 13,  201, 95,  96,  53,  194, 233, ...
                        7,   225, 140, 36,  103, 30,  69,  142, 8,   99,  37,  240, 21,  10, ...
                        23,  190, 6,   148, 247, 120, 234, 75,  0,   26,  197, 62,  94,  252, ...
                        219, 203, 117, 35,  11,  32,  57,  177, 33,  88,  237, 149, 56,  87, ...
                        174, 20,  125, 136, 171, 168, 68,  175, 74,  165, 71,  134, 139, 48, ...
                        27,  166, 77,  146, 158, 231, 83,  111, 229, 122, 60,  211, 133, 230, ...
                        220, 105, 92,  41,  55,  46,  245, 40,  244, 102, 143, 54,  65,  25, ...
                        63,  161, 1,   216, 80,  73,  209, 76,  132, 187, 208, 89,  18,  169, ...
                        200, 196, 135, 130, 116, 188, 159, 86,  164, 100, 109, 198, 173, 186, ...
                        3,   64,  52,  217, 226, 250, 124, 123, 5,   202, 38,  147, 118, 126, ...
                        255, 82,  85,  212, 207, 206, 59,  227, 47,  16,  58,  17,  182, 189, ...
                        28,  42,  223, 183, 170, 213, 119, 248, 152, 2,   44,  154, 163, 70, ...
                        221, 153, 101, 155, 167, 43,  172, 9,   129, 22,  39,  253, 19,  98, ...
                        108, 110, 79,  113, 224, 232, 178, 185, 112, 104, 218, 246, 97,  228, ...
                        251, 34,  242, 193, 238, 210, 144, 12,  191, 179, 162, 241, 81,  51, ...
                        145, 235, 249, 14,  239, 107, 49,  192, 214, 31,  181, 199, 106, 157, ...
                        184, 84,  204, 176, 115, 121, 50,  45,  127, 4,   150, 254, 138, 236, ...
                        205, 93,  222, 114, 67,  29,  24,  72,  243, 141, 128, 195, 78,  66, ...
                        215, 61,  156, 180];
    elseif type == 'r'
        PerlinNoise = 0:1:255;
        PerlinNoise=shuffle(PerlinNoise);
    end
    
    PerlinNoise = [PerlinNoise PerlinNoise];
end

 function v=shuffle(v)
     v=v(randperm(length(v)));
 end
 

