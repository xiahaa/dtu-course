clc;close all;clear all;

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

addpath ../utils;

data_dir = '../data/EX_3_data/';

skip = [1 ];

if isempty(find(skip==1,1))
    %% exercise 1
    imnoisy = (imread(strcat(data_dir,'noisy_circles.png')));
    imreal = (imread(strcat(data_dir,'noise_free_circles.png')));

    f = [1;2;3];% label
    miuf = [70;130;190];% mean intensities
    alpha = 0.0005;% weight for likelihhod
    beta = 0.8;% weight for prior

    %% Q1
    configurationTruth = imreal;
    configurationTruth(configurationTruth == miuf(f(1))) = 1;
    configurationTruth(configurationTruth == miuf(f(2))) = 2;
    configurationTruth(configurationTruth == miuf(f(3))) = 3;

    buildHistogram(imnoisy,configurationTruth,f);% histogram for real configuration
    Etrue = calcLikelihood(miuf,imnoisy,configurationTruth,alpha) + calcSmoothnessPrior(configurationTruth, beta);% energy

    %% Q2
    % test = [1 0 1 0 0;0 1 1 1 2;1 0 0 3 5;1 0 1 5 1;8 2 3 1 1;];
    % V2 = calcSmoothnessPrior(im1, 1);
    % V21 = burteForceSol(im1, 1);

    figure;
    subplot(3,2,1);imshow(imreal); title('Raw','FontName','Arial','FontSize',20);
    subplot(3,2,2);imshow(imnoisy); title('Noisy','FontName','Arial','FontSize',20);

    % configuration 1, simple threshold
    threshold1 = 100;
    threshold2 = 160;
    configuration1 = imnoisy;
    configuration1(imnoisy<100) = f(1);
    configuration1(imnoisy>=100 & imnoisy <= 160) = f(2);
    configuration1((imnoisy>=160)) = f(3);

    imseg1 = segmentation(imnoisy, configuration1, f, miuf);
    L1 = calcLikelihood(miuf,imnoisy,configuration1,alpha);
    P1 = calcSmoothnessPrior(configuration1, beta);
    E1 = L1 + P1;
    subplot(3,2,3);imshow(imseg1); title('Threshold','FontName','Arial','FontSize',20);

    % configuration 2, median filter on configuration 1
    configuration2 = medfilt2(configuration1,[3,3],'symmetric');
    imseg2 = segmentation(imnoisy, configuration2, f, miuf);
    L2 = calcLikelihood(miuf,imnoisy,configuration2,alpha);
    P2 = calcSmoothnessPrior(configuration2, beta);
    E2 = L2 + P2;
    subplot(3,2,4);imshow(imseg2); title('3x3 Median Filter','FontName','Arial','FontSize',20);

    % configuration 3, firstly gaussian smoothing and then take the
    % threshold
    imsmmoth = imgaussfilt(imnoisy,1,'FilterSize',3,'Padding','replicate');
    configuration3(imsmmoth<100) = f(1);
    configuration3(imsmmoth>=100 & imsmmoth <= 160) = f(2);
    configuration3((imsmmoth>=160)) = f(3);

    imseg3 = segmentation(imnoisy, configuration3, f, miuf);
    L3 = calcLikelihood(miuf,imnoisy,configuration3,alpha);
    P3 = calcSmoothnessPrior(configuration3, beta);
    E3 = L3 + P3
    subplot(3,2,5);imshow(imseg3); title('Gaussian','FontName','Arial','FontSize',20);

    % configuration 4, apply morphological operation (open) on
    % configuration1
    se = strel('disk',1);
    configuration4 = imopen(configuration1,se);
    imseg4 = segmentation(imnoisy, configuration4, f, miuf);
    L4 = calcLikelihood(miuf,imnoisy,configuration4,alpha);
    P4 = calcSmoothnessPrior(configuration4, beta);
    E4 = L4 + P4;
    subplot(3,2,6);imshow(imseg4); title('Morpho','FontName','Arial','FontSize',20);
    
    % display results
    figure
    subplot(2,2,1);
    buildHistogram(imnoisy,configuration1,f); title('Threshold','FontName','Arial','FontSize',20);
    subplot(2,2,2);
    buildHistogram(imnoisy,configuration2,f); title('3x3 Median Filter','FontName','Arial','FontSize',20);
    subplot(2,2,3);
    buildHistogram(imnoisy,configuration3,f); title('Gaussian filter','FontName','Arial','FontSize',20);
    subplot(2,2,4);
    buildHistogram(imnoisy,configuration4,f); title('Morphological: Close','FontName','Arial','FontSize',20);

    % plot likelihood cost
    figure
    bar([L1 L2 L3 L4]);
    xticklabels({'CF1','CF2','CF3','CF4'});
    xtickangle(45);
    xlim=get(gca,'xlim');
    hold on
    title('Likelihood');
    set(gca,'FontName','Arial','FontSize',20);

    % plot smoothness prior cost
    figure
    bar([P1 P2 P3 P4]);
    xticklabels({'CF1','CF2','CF3','CF4'});
    xtickangle(45);
    xlim=get(gca,'xlim');
    hold on
    title('Prior');
    set(gca,'FontName','Arial','FontSize',20);

    % plot posterior cost
    figure
    bar([E1 E2 E3 E4]);
    xticklabels({'CF1','CF2','CF3','CF4'});
    xtickangle(45);
    set(gca,'FontName','Arial','FontSize',20);
    xlim=get(gca,'xlim');
    hold on
    h1 = plot(xlim,[Etrue Etrue],'r--');
    legend([h1],{'Truth'});
    title('Posterior');

    %% Q3-ICM
    iter = 1;% current iteration
    maxiter = 1e4;% maximum iteration count
    oldCost = 1e6;% cost of previous iteration

    %configuration = configuration1;% start with cfg1
    configuration = randi(3,size(imnoisy));% start with random configuration

    % create a checkerboard with 1 pixel width for partial parallel
    % optimization
    mask = checkerboard(1,round(size(imnoisy,1)*0.5),round(size(imnoisy,2)*0.5));
    mask(size(imnoisy,1)+1:end,:) = [];
    mask(:,size(imnoisy,2)+1:end) = [];
    mask(mask > 0.5) = 1; mask = mask == 1;

    T = 1e4;% temprature if using simulated anneling

    figure
    
    optType = 3;
    optHandler = {@ICM,@GibbsSampling,@SimulatedAnneling};
    optName = {'ICM','Gibbs','SA'};
        
    % optmization loop
    while iter < maxiter
        fhandle = optHandler{optType};
        
        if strcmp(optName{optType},'SA')
            configuration = fhandle(imnoisy, configuration, f, miuf, alpha, beta, mask,T);
            T = T * 0.8;
        else
            configuration = fhandle(imnoisy, configuration, f, miuf, alpha, beta, mask);
        end

        % new cost
        Lnew = calcLikelihood(miuf,imnoisy,configuration,alpha);
        Pnew = calcSmoothnessPrior(configuration, beta);
        newCost = Lnew + Pnew;
    
        if abs(newCost-oldCost) < 1e-6
            break;
        end
        
        oldCost = newCost;
        costs(iter) = newCost;
        iter = iter + 1;
        % show variation
        imres = segmentation(imnoisy, configuration, f, miuf);
        imshow(imres);pause(0.1);
    end
    
    % cost cost variation with iterations
    figure
    plot(costs,'r-o','LineWidth',2);grid on;
    xlabel('Iteration');
    ylabel('Cost')
    title(optName{optType})
    set(gca,'FontName','Arial','FontSize',20);

    % plot posterior cost
    figure
    bar([E1 E2 E3 E4 newCost]);
    xticklabels({'CF1','CF2','CF3','CF4',optName{optType}});
    xtickangle(45);
    set(gca,'FontName','Arial','FontSize',20);
    xlim=get(gca,'xlim');
    hold on
    h1 = plot(xlim,[Etrue Etrue],'r--');
    legend([h1],{'Truth'});
    title('Posterior');
    
    %% use multilabel cut
    localPotential = calclocalPotentials(imnoisy, [], f, miuf, alpha, beta);
    U = reshape(localPotential,size(localPotential,1)*size(localPotential,2),size(localPotential,3));
    dim = [size(localPotential,1),size(localPotential,2)];
    max_iter = 1e6;
    [S,iter] = multilabel_MRF(U,dim,beta,max_iter);
    
    % show variation
    imres = segmentation(imnoisy, S, f, miuf);
    figure;imshow(imres);
    
    Lnew = calcLikelihood(miuf,imnoisy,S,alpha);
    Pnew = calcSmoothnessPrior(S, beta);
    newCostCut = Lnew + Pnew;
    
    % plot posterior cost
    figure
    bar([E1 E2 E3 E4 newCost newCostCut]);
    xticklabels({'CF1','CF2','CF3','CF4',optName{optType},'Cut'});
    xtickangle(45);
    set(gca,'FontName','Arial','FontSize',20);
    xlim=get(gca,'xlim');
    hold on
    h1 = plot(xlim,[Etrue Etrue],'r--','LineWidth',2);
    legend([h1],{'Truth'});
    title('Posterior');
    
end

if isempty(find(skip==2,1))
    %% exercise 1
    im1 = (imread(strcat(data_dir,'V12_10X_x502.png')));% 
%     im1 = (imread(strcat(data_dir,'V8_10X_x502.png')));%
    im1 = im2double(im1);
%     im1 = imresize(im1,0.3);
    
    figure;imshow(im1);
    
    figure;
    buildHistogram(im1,[],[]);
    
    f = [1,2];
    miuf = [0.4,0.7];
    
    addpath(genpath('./'));% addpath for graphcut 
    
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
    beta = 0.1; % weight of the prior term
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
    imres = segmentation(im1, configuration, f, miuf);
    figure;imshow(imres);
    
    figure
    buildHistogram(im1,configuration,f);
    
    %% use multilabel cut
    fmulti = [1 2 3];
    miufulti = [0.4,0.45,0.7];
    alpha = 1;
    beta = 0.015;
    
    localPotential = calclocalPotentials(im1, [], fmulti, miufulti, alpha, beta);
    U = reshape(localPotential,size(localPotential,1)*size(localPotential,2),size(localPotential,3));
    dim = [size(localPotential,1),size(localPotential,2)];
    max_iter = 1e6;
    [S,iter] = multilabel_MRF(U,dim,beta,max_iter);
    
    % show variation
    ccmap = jet(numel(fmulti));
    imseg1 = zeros(size(localPotential,1),size(localPotential,2));
    imseg2 = zeros(size(localPotential,1),size(localPotential,2));
    imseg3 = zeros(size(localPotential,1),size(localPotential,2));
    for i = 1:length(fmulti)
        imseg1(S == fmulti(i)) = ccmap(i,1);
        imseg2(S == fmulti(i)) = ccmap(i,2);
        imseg3(S == fmulti(i)) = ccmap(i,3);
    end
    imseg = cat(3,imseg1,imseg2,imseg3);
    figure;imshow(imseg);
    
    figure
    buildHistogram(im1,S,fmulti);
    
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
