clc;close all;clear all;

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

% add path to the data directory
addpath ../data/EX_2_data;
addpath ../utils/

skip = [1 2 ];
drawGif = 0;

if isempty(find(skip == 1,1))
    %% ex1
    N = 2;
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
    plot(p(1,:),p(2,:),'r*','MarkerSize',10);
    hold on;
    plot(q(1,:),q(2,:),'bo','MarkerSize',10);
    title('Inputs','FontSize',fontsize);
    legend({'p','q'},'FontSize',fontsize);
    xlabel('x: (m)','FontSize',fontsize);
    ylabel('y: (m)','FontSize',fontsize);
    axis equal;grid on;

    [R1,t1,s1] = lsq_allign(p,q);
    qr = R1'*(q./s1-repmat(t1,1,N));
    subplot(1,2,2);
    plot(p(1,:),p(2,:),'r*','MarkerSize',10);
    hold on;
    plot(qr(1,:),qr(2,:),'bo','MarkerSize',10);
    title('Recover','FontSize',fontsize);
    legend({'p','q'},'FontSize',fontsize);
    xlabel('x: (m)','FontSize',fontsize);
    ylabel('y: (m)','FontSize',fontsize);
    axis equal;grid on;

    set(gca,'fontname','arial')  % Set it to Arial
    
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
    plot(p(1,:),p(2,:),'r*','MarkerSize',10);
    hold on;
    plot(q(1,:),q(2,:),'bo','MarkerSize',10);
    title('Inputs','FontSize',fontsize);
    legend({'p','q'},'FontSize',fontsize);
    xlabel('x: (m)','FontSize',fontsize);
    ylabel('y: (m)','FontSize',fontsize);
    axis equal;grid on;

    [R1,t1,s1] = lsq_allign(p,q);
    qr = R1'*(q./s1-repmat(t1,1,N));
    subplot(1,2,2);
    plot(p(1,:),p(2,:),'r*','MarkerSize',10);
    hold on;
    plot(qr(1,:),qr(2,:),'bo','MarkerSize',10);
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
    shift = size(im1,2);
    cim(1:size(im1,1),1:size(im1,2)) = im1;
    cim(1:size(im2,1),shift+1:shift+size(im2,2)) = im2;
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
    ransac.minimum_sample = 2;
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
        
        if isempty(R)
            iter = iter + 1;
            continue;
        end
        
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

    inlier_maches = matches(ransac.inliers,:);
    
    figure; imshow(cim); hold on;
    h1 = vl_plotframe(feature1(:,matches(:,1))) ;
    set(h1,'color','g','linewidth',1) ;
    for i = 1:size(inlier_maches,1)
        xx = [feature1(1,inlier_maches(i,1)) feature2(1,inlier_maches(i,2))+shift];
        yy = [feature1(2,inlier_maches(i,1)) feature2(2,inlier_maches(i,2))];
        plot(xx,yy,'r-');
    end
    fontsize = 15;
    title('Refined Matching','FontSize',fontsize);
    set(gca,'fontname','arial')  % Set it to Arial
    
    figure
    imshow(im2); hold on;
    h2 = vl_plotframe(feature2(:,inlier_maches(:,2)));
    set(h2,'color','b','linewidth',1) ;
    plot(q(2,ransac.inliers),q(1,ransac.inliers),'co');
    plot(qt(2,ransac.inliers),qt(1,ransac.inliers),'m*');
    title('Alignment of 2D features','FontSize',fontsize);
    set(gca,'fontname','arial')  % Set it to Arial
end

if isempty(find(skip == 3,1))
    vl_setup;

    im1id = 1;
    im2id = 2;
    
    im1 = imread(imgnames{im1id});
    im2 = imread(imgnames{im2id});
    
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
    shift = size(im1,2);
    cim(1:size(im1,1),1:size(im1,2)) = im1;
    cim(1:size(im2,1),shift+1:shift+size(im2,2)) = im2;
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
    ransac.minimum_sample = 2;
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
        
        if isempty(R)
            iter = iter + 1;
            continue;
        end
        
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
    
    inlier_maches = matches(ransac.inliers,:);

    figure
    imshow(im2); hold on;
    h2 = vl_plotframe(feature2(:,inlier_maches(:,2)));
    set(h2,'color','b','linewidth',1) ;
    plot(q(2,ransac.inliers),q(1,ransac.inliers),'co');
    plot(qt(2,ransac.inliers),qt(1,ransac.inliers),'m*');
    
    smooth_scales = [4 3 1 4 4 4];

    t1 = smooth_scales(im1id);%1 4 2 3 3 1 4 4, 5 4,6 4
    sigma = 3;
    g1 = gassian_fast(t1, sigma);
    im1blur = imfilter(im1, g1, 'replicate');
    im1blur = imfilter(im1blur, g1', 'replicate');
    
    t2 = smooth_scales(im2id);%1 4 2 3 3 1 4 4, 5 4,6 4
    g2 = gassian_fast(t2, sigma);
    im2blur = imfilter(im2, g2, 'replicate');
    im2blur = imfilter(im2blur, g2', 'replicate');
    
    [local_maxima1, rs1] = multi_scale_detect_blob(im1, im1blur);
    [local_maxima2, rs2] = multi_scale_detect_blob(im2, im2blur);
    local_maxima1 = local_maxima1';
    local_maxima2 = local_maxima2';
    
    figure
    imshow(cim); hold on;
    w = linspace(0,2*pi,100);
    for i = 1:size(local_maxima1,2)
        cx = local_maxima1(2,i);
        cy = local_maxima1(1,i);
        xx = cx + rs1(i).*cos(w);
        yy = cy + rs1(i).*sin(w);
        plot(xx,yy,'b-');
    end
    for i = 1:size(local_maxima2,2)
        cx = local_maxima2(2,i);
        cy = local_maxima2(1,i);
        xx = cx + rs2(i).*cos(w) + shift;
        yy = cy + rs2(i).*sin(w);
        plot(xx,yy,'r-');
    end
    
    local_maxima2to1 = R'*(local_maxima2./s-repmat(t,1,size(local_maxima2,2)));
    
    % 2D to 1D
    [M,N] = size(im1);
    indices_maxima2to1 = local_maxima2to1(1,:) + local_maxima2to1(2,:).*M;
    indices_maxima1 = local_maxima1(1,:) + local_maxima1(2,:).*M;
    
    n1 = size(indices_maxima2to1,2);
    n2 = size(indices_maxima1,2);
    indices_maxima2to1 = repmat(indices_maxima2to1,n2,1);
    indices_maxima1 = repmat(indices_maxima1',1,n1);
    diff_indices = indices_maxima2to1 - indices_maxima1;
    [mindiff,minid] = min(abs(diff_indices),[],1);
    matched = mindiff < M+1;
    matches12 = zeros(n1,2);
    matches12(matched,1) = minid(matched)';
    tmp = 1:n1;
    matches12(matched,2) = tmp(matched)';
    
%     s.*(R*local_maxima1+repmat(t,1,size(local_maxima1,2)));
%     matches12 = zeros(size(local_maxima2to1,2),2);
%     for i = 1:size(local_maxima2to1,2)
%         diff12 = repmat(local_maxima2to1(:,i),1,size(local_maxima1,2))-local_maxima1;
%         err = diag(diff12'*diff12);
%         [minval,minid] = min(err);
%         if minval > 4
%             continue;
%         end
%         matches12(i,:) = [minid, i];
%     end
    local_maxima2to1 = local_maxima2to1(:,matches12(:,2) ~= 0);

    figure
    imshow(im1); hold on;
    w = linspace(0,2*pi,100);
    rs2 = rs2./s;
    for i = 1:size(local_maxima2to1,2)
        cx = local_maxima2to1(2,i);
        cy = local_maxima2to1(1,i);
        xx = cx + rs2(i).*cos(w);
        yy = cy + rs2(i).*sin(w);
        plot(xx,yy,'r-');
    end
    for i = 1:size(local_maxima1,2)
        cx = local_maxima1(2,i);
        cy = local_maxima1(1,i);
        xx = cx + rs1(i).*cos(w);
        yy = cy + rs1(i).*sin(w);
        plot(xx,yy,'b-');
    end
    
    % some statistics
    disp(['mean diameters 1: ', num2str(mean(rs1))]);
    disp(['mean diameters 2: ', num2str(mean(rs2))]);
    
    disp(['std of diameters 1: ', num2str(std(rs1))]);
    disp(['std of diameters 2: ', num2str(std(rs2))]);
    
    validid = matches12(:,2) ~= 0;
    diff_of_matches = abs(rs1(matches12(validid,1)) - rs2(matches12(validid,2)));
    figure;hist(diff_of_matches);
end



function [local_maxima, rs] = multi_scale_detect_blob(Igray, Iblur)
    half_win_size = 1;
    [local_maxima, ~] = find_maxima(Iblur, half_win_size);
    K = 6;
    t0 = 2^(1/5);
    rs = zeros(size(local_maxima,1),1);
    [LLNs, radius] = create_scale_normalized_LoG(Igray, t0, K);
    for i = 1:size(local_maxima,1)
        val = LLNs(local_maxima(i,1),local_maxima(i,2),:);
        [~,minid] = min((val));
        rs(i) = radius(minid);
    end
end
    
%% separate functions
function [R,t,s] = lsq_allign(p,q)
    qm = mean(q,2);
    pm = mean(p,2);
    qd = q - repmat(qm,1,size(q,2));
    pd = p - repmat(pm,1,size(p,2));
    
    % two match one case
    if all(abs(qd(:))) < 1e-6 || all(abs(pd(:))) < 1e-6
        R=[];t=[];s=[];
        return;
    end
    
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
        R = U*diag([1,-1])*V';
    end
    t = qm - R*pm;
end