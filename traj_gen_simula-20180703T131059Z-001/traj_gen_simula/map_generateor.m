function [output, pcl] = map_generateor(sizeX, sizeY, sizeZ, seed, scale)
    clc;clear all;
    
    resolution = 0.25;%% 1 grid = 0.25m
    scale = 1/resolution;
    sizeX = 100;
    sizeY = 100;
    sizeZ = 5;
    
%     sizeX = sizeX * resolution;%from grid number to real size
%     sizeY = sizeY * resolution;
%     sizeZ = sizeZ * resolution;
    
    %% obstacles, w_l = minimium width of obstacles in real metrics
    ObsNum = 30;
    w_l = 1;
    w_h = 2;
    
    map3D = robotics.OccupancyMap3D(scale);
    pose = [ 0 0 0 1 0 0 0];
    maxRange = sizeX * 0.6;
    
%     perlin3D();
    pclXYZ = randomMapGen(sizeX, sizeY, sizeZ, scale, ObsNum, w_l, w_h);
%     pclXYZ = perlin3D(sizeX, sizeY, sizeZ, scale);
    
    insertPointCloud(map3D,pose,pclXYZ,maxRange);
    
    show(map3D);
    grid on;
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
    complexity = 0.142857;
    fill = 0.38;
    fractal = 1;
    attenuation = 0.5;
    
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
                    tnoise = tnoise + ta * calcNoise(dfv * i * complexity, dfv * j * complexity, dfv * k * complexity, noise);
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
                    tnoise = tnoise + ta * calcNoise(dfv * i * complexity, dfv * j * complexity, dfv * k * complexity, noise);
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
 
 function res = calcNoise(x,y,z,p)
    X = bitand(floor(x),255,'uint8');
    Y = bitand(floor(y),255,'uint8');
    Z = bitand(floor(z),255,'uint8');
    
    x = x-floor(x);
    y = y-floor(y);
    z = z-floor(z);

    % Compute fade curves for each of x, y, z
    u = fade(x);
    v = fade(y);
    w = fade(z);
    
    % Hash coordinates of the 8 cube corners
    A  = p(X+1) + Y;
    AA = p(A+1) + Z;
    AB = p(A + 2) + Z;
    B  = p(X + 2) + Y;
    BA = p(B) + Z;
    BB = p(B + 1) + Z;
    
    % Add blended results from 8 corners of cube
    res = lerp(w, lerp(v, lerp(u, grad(p(AA+1), x, y, z), grad(p(BA+1), x - 1, y, z)), ...
              lerp(u, grad(p(AB+1), x, y - 1, z), grad(p(BB+1), x - 1, y - 1, z))), ...
                lerp(v, lerp(u, grad(p(AA+2), x, y, z - 1), grad(p(BA + 2), x - 1, y, z - 1)), ...
                    lerp(u, grad(p(AB + 2), x, y - 1, z - 1), grad(p(BB + 2), x - 1, y - 1, z - 1))));
    res = (res + 1.0) / 2.0;
 end

function n = fade(t)
    n = t * t * t * (t * (t * 6 - 15) + 10);
end

function n = lerp(t,a,b)
    n = a + t * (b - a);
end

function n = grad(hash,x,y,z)
    h = bitand(hash,15,'int32');%    hash & 15;
    if h < 8
        u = x;
    else
        u = y;
    end
    if h < 4
        v = y;
    else
        if h == 12 || h == 14
            v = x;
        else
            v = z;
        end
    end
    if bitand(h,1,'int32') == 0
        if bitand(h,2,'int32') == 0
            n = u + v;
        else
            n = u - v;
        end
    else
        if bitand(h,2,'int32') == 0
            n = -u + v;
        else
            n = -u - v;
        end
    end
end