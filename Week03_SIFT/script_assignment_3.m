clc;close all;clear all;

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

% add path to the data directory
addpath ../data/EX_2_data;
addpath ../utils/

skip = [1 ];
drawGif = 0;

if isempty(find(skip == 1,1))
    %% ex1
    N = 100;
    dim = 2;
    p = rand(dim,N);
    theta = rand(1)*(2*pi)-pi;% -pi to pi
    R = [cos(theta) sin(theta);-sin(theta) cos(theta)];
    t = rand(2,1);
    s = rand(1)*0.9 + 0.1;
    q = (R*p + repmat(t,1,N));
    q = s.*q;
    figure
    subplot(1,2,1);fontsize = 15;
    plot(p(1,:),p(2,:),'r*','MarkerSize',5);
    hold on;
    plot(q(1,:),q(2,:),'bo','MarkerSize',5);
    title('Inputs','FontSize',fontsize);
    legend({'p','q'},'FontSize',fontsize);
    xlabel('x: (m)','FontSize',fontsize);
    ylabel('y: (m)','FontSize',fontsize);
    axis equal;grid on;

    [R1,t1,s1] = lsq_allign(p,q);
    qr = R1'*(q./s1-repmat(t1,1,N));
    subplot(1,2,2);
    plot(p(1,:),p(2,:),'r*','MarkerSize',5);
    hold on;
    plot(qr(1,:),qr(2,:),'bo','MarkerSize',5);
    title('Recover','FontSize',fontsize);
    legend({'p','q'},'FontSize',fontsize);
    xlabel('x: (m)','FontSize',fontsize);
    ylabel('y: (m)','FontSize',fontsize);
    axis equal;grid on;

    set(gca,'fontname','arial')  % Set it to Arial
end

imgnames = {'CT_lab_high_res.png', ...
             'CT_lab_med_res.png', ...
             'CT_lab_low_res.png', ...
             'CT_synchrotron.png', ...
             'Optical.png', ...
             'SEM.png'};

if isempty(find(skip == 2,1))
    vl_setup;

    im1 = imread(imgnames{1});
    im2 = imread(imgnames{2});
    
    % rgb2gray
    if size(im1,3) == 3
        im1 = rgb2gray(im1);
    end
    if size(im2,3) == 3
        im2 = rgb2gray(im2);
    end
    
    % to single for vlfeat
    im1single = single(im1);
    im2single = single(im2);
    
    im1 = im2double(im1);
    im2 = im2double(im2);
    
    % feature detection and description
    [feature1,desp1] = vl_sift(im1single);
    [feature2,desp2] = vl_sift(im2single);
    
    desp1 = single(desp1);desp2 = single(desp2);
    % start matching using kd-tree
    kdtree = vl_kdtreebuild(desp2);
    threshold = 1.6;
    % pre-allocate
    matches = zeros(size(feature1,2),2);
    k = 1;
    for i = 1:size(feature1,2)
        [index, distance] = vl_kdtreequery(kdtree,desp2,desp1(:,i),'NumNeighbors',2, 'MaxComparisons', 15);
        if threshold*distance(1) < distance(2)
            matches(k,:) = [i, index(1)];
            k = k + 1;
        end
    end
    % delete unused memory
    matches(k:end,:) = [];
    
    % show
    cwidth = size(im1,2) + size(im2,2);
    cheight = max([size(im1,1),size(im2,1)]);
    cim = zeros(cheight,cwidth);
    shift = size(im1,1);
    cim(1:size(im1,1),1:size(im1,2)) = im1;
    cim(1:size(im2,2),shift+1:shift+size(im2,2)) = im2;
    figure; imshow(cim); hold on;
    h1 = vl_plotframe(feature1(:,matches(:,1))) ;
    set(h1,'color','g','linewidth',1) ;
    for i = 1:size(matches,1)
        xx = [feature1(1,matches(i,1)) feature2(1,matches(i,2))+shift];
        yy = [feature1(2,matches(i,1)) feature2(2,matches(i,2))];
        plot(xx,yy,'r-');
    end
    fontsize = 15;
    title('Raw Matching','FontSize',fontsize);
    set(gca,'fontname','arial')  % Set it to Arial
    
    ransac.maxiter = 1e6;
    ransac.threshold = 0.25;% within 0.5 pixel, then it is a inlier
    ransac.prob = 0.99;
    ransac.minimum_sample = 3;
    ransac.bestcost = 0;
    ransac.inliers = [];
    iter = 1;
    
    % data preparation
    p = [feature1(2,matches(:,1));feature1(1,matches(:,1))];
    q = [feature2(2,matches(:,2));feature2(1,matches(:,2))];
    N = size(matches,1);
    while iter<ransac.maxiter
        samples = randperm(N,ransac.minimum_sample);
        ps = p(:,samples);
        qs = q(:,samples);
        [R,t,s] = lsq_allign(ps,qs);
        
        % transform
        qt = s.*(R*p+repmat(t,1,N));
        err = qt - q;
        err = diag(err'*err)';
        % inliers
        inlier_id = err < ransac.threshold;
        inlier_cnt = sum(inlier_id);
        % update
        if ransac.bestcost < inlier_cnt
            ransac.inliers = inlier_id;
            ransac.bestcost = inlier_cnt;
            inlier_percentage = inlier_cnt / N;
            ransac.maxiter = log(0.01)/log(1-inlier_percentage^ransac.minimum_sample);
        end
        iter = iter + 1;
    end
    %% refine with inliers
    ps = p(:,ransac.inliers);
    qs = q(:,ransac.inliers);
    [R,t,s] = lsq_allign(ps,qs);
    qt = s.*(R*p+repmat(t,1,N));
    
    figure; imshow(cim); hold on;
    h1 = vl_plotframe(feature1(:,matches(:,1))) ;
    set(h1,'color','g','linewidth',1) ;
    inlier_maches = matches(ransac.inliers,:);
    for i = 1:size(inlier_maches,1)
        xx = [feature1(1,inlier_maches(i,1)) feature2(1,inlier_maches(i,2))+shift];
        yy = [feature1(2,inlier_maches(i,1)) feature2(2,inlier_maches(i,2))];
        plot(xx,yy,'r-');
    end
    fontsize = 15;
    title('Raw Matching','FontSize',fontsize);
    set(gca,'fontname','arial')  % Set it to Arial
    
    figure
    imshow(im2); hold on;
    h2 = vl_plotframe(feature2(:,inlier_maches(:,2)));
    set(h2,'color','b','linewidth',1) ;
    plot(q(2,ransac.inliers),q(1,ransac.inliers),'bo');
    plot(qt(2,ransac.inliers),qt(1,ransac.inliers),'m*');
end



%% separate functions
function [R,t,s] = lsq_allign(p,q)
    qm = mean(q,2);
    pm = mean(p,2);
    qd = q - repmat(qm,1,size(q,2));
    pd = p - repmat(pm,1,size(p,2));
    s = sqrt(sum(diag(qd'*qd))/sum(diag(pd'*pd)));
    qs = qd./s;
    qm = qm./s;
    C = zeros(2,2);
    for i = 1:size(p,2)
        C = C + qs(:,i)*pd(:,i)';
    end
    [U,~,V] = svd(C);
    R = U*V';
    if det(R)<0
        R = U*diag(1,-1)*V';
    end
    t = qm - R*pm;
end