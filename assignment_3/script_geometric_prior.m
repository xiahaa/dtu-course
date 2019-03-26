clc;close all;clear all;

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

addpath ../utils;
addpath(genpath('./'));% addpath for graphcut 

data_dir = '../data/EX_5_data/';

skip = [1 2 ];

%% exercise 1
info = imfinfo(strcat(data_dir,'nerves_part.tiff'));
numOfImages = size(info,1);
    
    if isempty(find(skip==1,1))
        im1 = imread(strcat(data_dir,'nerves_part.tiff'),1);% 

        f = [1,2];
        miuf = [0.15,0.5];
    
        gen3D = 0;
    
        if gen3D == 0
            %% Q1, slice by slice
            truncatedNum = 700;
            confs = zeros(size(im1,1),size(im1,2),truncatedNum);
            for k = 1:1:numOfImages
                disp(k);
                im1 = imread(strcat(data_dir,'nerves_part.tiff'),k);% 
                figure(1);subplot(1,2,1);imshow(im1);
                im1 = im2double(im1);

        %         figure;
        %         buildHistogram(im1,[],[]);

                d = im1(:); % intensity (data)
                mu = miuf; % means of two classes

                %% establishing likelihood
                w_s = (d(:)-mu(1)).^2; % source weight
                w_t = (d(:)-mu(2)).^2; % sink weights
                N = numel(d); % number of graph nodes
                indices = (1:N)'; % an index for each person
                % terminal edge matrix
                E_terminal = [indices,[w_s,w_t]]; 

                %% establishing prior
                beta = 0.015; % weight of the prior term
                % internal edge matrix
                m = size(im1,1);n = size(im1,2);

                nn1 = zeros(n*(m-1),2);
                % row sweeping
                c1 = ((1:n)-1).*m; c1 = c1';
                for i = 1:m-1   
                    nn1((i-1)*n+1:i*n,:) = [c1+i c1+i+1];
                end
                % col sweeping
                c1 = 1:m; c1 = c1';
                nn2 = zeros(m*(n-1),2);
                for i = 1:n-1   
                    nn2((i-1)*m+1:i*m,:) = [c1+(i-1)*m c1+(i)*m];
                end

                % all neighbors
                nn = [nn1;nn2];
                E_internal = [nn(:,1),nn(:,2),beta*ones(size(nn,1),2)]; 

                [Scut,flow] = GraphCutMex(N,E_terminal,E_internal); % here it happens

                S = f(1).*ones(1,N);
                S(Scut) = f(2);
                configuration = reshape(S,m,n);

                % show variation
                ccmap = lines(numel(f));
                imseg1 = zeros(m,n);
                imseg2 = zeros(m,n);
                imseg3 = zeros(m,n);
                for i = 1:length(f)
                    imseg1(configuration == f(i)) = ccmap(i,1);
                    imseg2(configuration == f(i)) = ccmap(i,2);
                    imseg3(configuration == f(i)) = ccmap(i,3);
                end
                imseg = cat(3,imseg1,imseg2,imseg3);
                figure(1);subplot(1,2,2);imshow(imseg);

%                 F = getframe(gcf);
%                 imwrite(F.cdata,strcat('../data/Ex_5_data/output_mrf_1/',num2str(k,'%d'),'.png'));
                confs(:,:,k) = configuration;
                if k == truncatedNum
                    break;
                end
            end
            save('mrf_1.mat','confs');
        else
            if gen3D == 1
                load('mrf_1.mat');
                numOfImages = size(confs,3);
                V = zeros(size(im1,1),size(im1,2),numOfImages);
                [X,Y,Z] = meshgrid(1:size(im1,1),1:size(im1,2),1:numOfImages);
                for i = 1:numOfImages
                    imn = zeros(size(im1,1),size(im1,2));
                    imn(confs(:,:,i) == f(1)) = 1;
                    V(:,:,i) = (imn);
                end
            end
        end

        if gen3D == 1
            save('MRF.mat','V');

            pp = patch(isosurface(X,Y,Z,V,0));
            isonormals(X,Y,Z,V,pp)
            pp.FaceColor = [153/255,0/255,0/255];
            pp.EdgeColor = 'none';
            daspect([1 1 1]);
            view(3);
            axis tight;
            camlight;
            lighting gouraud
            title('Volumetric Image','Interpreter','latex','FontSize',15);
        end
    end
   
    %% Q2: full 3D segmentation
    if isempty(find(skip==2,1))
        f = [1,2];
        miuf = [0.2,0.5];
        mu = miuf; % means of two classes
        im1 = imread(strcat(data_dir,'nerves_part.tiff'),1);% 
        m = size(im1,1);n = size(im1,2);
        N = numel(im1); % number of graph nodes
        d = 100;
        
        startid = 0;
        shouldExit = 0;
        
        configurationfinal = zeros(m,n,numOfImages);
        
        doSegmentation = 0;
        
        if doSegmentation == 1
            while shouldExit == 0
                if (startid + d) > numOfImages
                    d = numOfImages - startid;
                    shouldExit = 1;
                end

                E_terminal = zeros(m*n*d, 1+numel(f));
                for k = 1:d
                    kk = k + startid;
                    disp(kk);
                    im1 = imread(strcat(data_dir,'nerves_part.tiff'),kk);% 
                    figure(1);imshow(im1);
                    im1 = im2double(im1);

                    %% establishing likelihood
                    w_s = (im1(:)-mu(1)).^2; % source weight
                    w_t = (im1(:)-mu(2)).^2; % sink weights

                    indices = (1:N)'+(k-1)*N; % an index for each person
                    % terminal edge matrix
                    E_terminal(indices,:) = [indices,[w_s,w_t]]; 
                end
                %% establishing prior
                beta = 0.015; % weight of the prior term
                dim = [m,n,d];

                indices = reshape(1:prod(dim),dim);
                edge_x = indices(1:end-1,:,:);
                edge_y = indices(:,1:end-1,:);
                edge_z = indices(:,:,1:end-1); % empty for 2D segmentation
                edge_n = zeros(numel(edge_x)+numel(edge_y)+numel(edge_z),4);
                edge_n(:,[1,2]) = [edge_x(:),edge_x(:)+1; edge_y(:),edge_y(:)+dim(1); ...
                    edge_z(:),edge_z(:)+dim(1)*dim(2)];

                edge_n(:,[3,4]) = beta;

                [Scut,flow] = GraphCutMex(m*n*d,E_terminal,edge_n); % here it happens
                S = f(1).*ones(1,m*n*d);
                S(Scut) = f(2);
                configurations = reshape(S,m,n,d);

                % show variation
                for j = 1:d
                    configuration = configurations(:,:,j);
                    ccmap = lines(numel(f));
                    imseg1 = zeros(m,n);
                    imseg2 = zeros(m,n);
                    imseg3 = zeros(m,n);
                    for i = 1:length(f)
                        imseg1(configuration == f(i)) = ccmap(i,1);
                        imseg2(configuration == f(i)) = ccmap(i,2);
                        imseg3(configuration == f(i)) = ccmap(i,3);
                    end
                    imseg = cat(3,imseg1,imseg2,imseg3);
                    imshow(imseg);pause(0.1);
                end

                configurationfinal(:,:,startid+1:startid+d) = configurations;
                startid = startid + d;
            end
            save('MRF3.mat','configurationfinal');
        else
            load('MRF3.mat');
            
            numOfImages = 200;
            V = zeros(size(im1,1),size(im1,2),numOfImages);
            [X,Y,Z] = meshgrid(1:size(im1,1),1:size(im1,2),1:numOfImages);
            
            for i = 1:min(size(configurationfinal,3),numOfImages)
                im1 = imread(strcat(data_dir,'nerves_part.tiff'),i);% 
                configuration = configurationfinal(:,:,i);
                % show variation
                ccmap = lines(numel(f));
                imseg1 = zeros(m,n);
                imseg2 = zeros(m,n);
                imseg3 = zeros(m,n);
                for j = 1:length(f)
                    imseg1(configuration == f(j)) = ccmap(j,1);
                    imseg2(configuration == f(j)) = ccmap(j,2);
                    imseg3(configuration == f(j)) = ccmap(j,3);
                end
                imseg = cat(3,imseg1,imseg2,imseg3);
                figure(1);subplot(1,2,1);imshow(im1);
                figure(1);subplot(1,2,2);imshow(imseg);
                
                imn = zeros(size(im1,1),size(im1,2));
                imn(configuration(:,:) == f(1)) = 1;
                V(:,:,i) = (imn);
            end
            
            pp = patch(isosurface(X,Y,Z,V,0));
            isonormals(X,Y,Z,V,pp)
            pp.FaceColor = lines(1);
            pp.EdgeColor = 'none';
            daspect([1 1 1]);
            view(3);
            axis tight;
            camlight;
            lighting gouraud
            title('Volumetric Image','Interpreter','latex','FontSize',15);
        end
    end
    
    %% Q3: deformable detection on segmented image
    if isempty(find(skip==3,1))
        f = [1,2];
        miuf = [0.2,0.5];
        mu = miuf; % means of two classes
        im1 = imread(strcat(data_dir,'nerves_part.tiff'),1);% 
        m = size(im1,1);n = size(im1,2);
        
        load('MRF3.mat');
        load('nerves.mat');
        
        % parameters
        Num = 50;
        stepSize = 30;
        a = 0.2;
        b = 0.8;
        % regularization
        Bint = regularization(a, b, Num);
        
        numNerves = size(ncc,1);
        
        % curve initialization, for each initialize a deformable curve
        alpha = linspace(0,2*pi, Num);
        for i = 1:numNerves
            r = 30;
            nerves(i).curve(1,:) = ncc(i,1) + r*cos(alpha);
            nerves(i).curve(2,:) = ncc(i,2) + r*sin(alpha);
        end
        figure(1);
        figure(2);
        se = strel('disk',40);
        ccmap = lines(numNerves);
        
        numOfImages = min(numOfImages, 500);
        
        doDeformable = 1;
        
        if doDeformable == 1
            for k = 1:1:numOfImages
                disp(k);
                im1 = imread(strcat(data_dir,'nerves_part.tiff'),k);%
                im1 = im2double(im1);
                conf = configurationfinal(:,:,k);

                mask = conf == f(1);

    %             ccmap = lines(numel(f));
    %             imseg1 = zeros(m,n);
    %             imseg2 = zeros(m,n);
    %             imseg3 = zeros(m,n);
    %             for i = 1:length(f)
    %                 imseg1(conf == f(i)) = ccmap(i,1);
    %                 imseg2(conf == f(i)) = ccmap(i,2);
    %                 imseg3(conf == f(i)) = ccmap(i,3);
    %             end
    %             imseg = cat(3,imseg1,imseg2,imseg3);
    %             figure(1);imshow(imseg);

                im2 = im1;
                im2(mask) = miuf(f(1));
                im2(~mask) = miuf(f(2));

                figure(2);imshow(im1); hold on;

                for i = 1:numNerves
                    plot(nerves(i).curve(2,[1:end,1]),nerves(i).curve(1,[1:end,1]),'-','Color',ccmap(i,:),'LineWidth',2);
                    nerves(i).curve = defromableCurve(im1, nerves(i).curve, Num, stepSize, Bint);
                end

    %             F = getframe;
    %             imwrite(F.cdata,strcat('../data/Ex_5_data/output_1/',num2str(k,'%d'),'.png'));

%                 save(strcat('../data/Ex_5_data/nerves/',num2str(k,'%d'),'.mat'), 'nerves');
            end
        else
            V = zeros(size(im1,1),size(im1,2),numOfImages);
            [X,Y,Z] = meshgrid(1:size(im1,1),1:size(im1,2),1:numOfImages);
            
            for k = 1:1:numOfImages
                disp(k);
                im1 = imread(strcat(data_dir,'nerves_part.tiff'),k);%
                im1 = im2double(im1);
                
                load(strcat('../data/Ex_5_data/nerves/',num2str(k,'%d'),'.mat'));
                
                figure(1);imshow(im1); hold on;

                imn = zeros(m,n);
                for i = 1:numNerves
                    mask = poly2mask(nerves(i).curve(2,:),nerves(i).curve(1,:),m,n);
                    boundary = findBoundary(mask);
                    se = strel('disk',1);
                    boundary = imdilate(boundary, se);
                    boundary = boundary == 1;
                    imn(boundary) = 1;
                end
                figure(2);imshow(imn);
                V(:,:,k) = (imn);
            end
            
            pp = patch(isosurface(X,Y,Z,V,0));
            isonormals(X,Y,Z,V,pp)
            pp.FaceColor = lines(1);
            pp.EdgeColor = 'none';
            daspect([1 1 1]);
            view(3);
            axis tight;
            camlight;
            lighting gouraud
            title('Volumetric Image','Interpreter','latex','FontSize',15);
        end
        
    end

function curve = defromableCurve(im, curve, Num, stepSize, Bint)
    m = size(im,1);
    n = size(im,2);

    mask = poly2mask(curve(2,:), curve(1,:), m, n);
    boundary = findBoundary(mask);
    se = strel('disk',3);
    boundary = imdilate(boundary, se);
    
    % find mean intensities inside and outside
    [cin, cout] = meanIntensity(im, curve, boundary);

    % compute the displacement along the normal direction
    displacement = computeDisplacement(im, curve, cin, cout);

    % find final
    displacement = displacement.*stepSize;
    %         displacement = findDisplacement(im, curve, displacement);

    % draw normals
%     quiver(curve(2,:),curve(1,:),displacement(2,:),displacement(1,:));

    % update
    curve = (Bint\(curve+displacement)')';

    % reinterpolation
    curve = reInterpolate(curve,Num);
%     curve = suppressSelfIntersection(curve, m, n);
    curve = constraintCurve(curve, m, n);
end

function buildHistogram(im1,seg,f)
    [h,wout] = hist(double(im1(:)),256);
    bar(wout,h);hold on;grid on;
    ccmap = lines(numel(f));
    for i = 1:numel(f)
        [h1,wout1] = hist(double(im1(seg==f(i))),256);  
%     [h2,wout2] = hist(double(im1(seg==2)),256);
%     [h3,wout3] = hist(double(im1(seg==3)),256);
        plot(wout1(h1~=0),h1(h1~=0),'-','LineWidth',2,'Color',ccmap(i,:));
    end
%     
%     plot(wout1(h1~=0),h1(h1~=0),'-','LineWidth',2,'Color',ccmap(1,:));
%     plot(wout2(h2~=0),h2(h2~=0),'-','LineWidth',2,'Color',ccmap(2,:));
%     plot(wout3(h3~=0),h3(h3~=0),'-','LineWidth',2,'Color',ccmap(3,:));
    title('Histogram','FontName','Arial','FontSize',20);
end
