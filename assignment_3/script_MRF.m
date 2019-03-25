clc;close all;clear all;

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

addpath ../utils;

data_dir = '../data/EX_3_data/';

imnoisy = (imread(strcat(data_dir,'noisy_circles.png')));
imreal = (imread(strcat(data_dir,'noise_free_circles.png')));

f = [1;2;3];
miuf = [70;130;190];
alpha = 0.0005;
beta = 1;

%% Q1
configurationTruth = imreal;
configurationTruth(configurationTruth == miuf(f(1))) = 1;
configurationTruth(configurationTruth == miuf(f(2))) = 2;
configurationTruth(configurationTruth == miuf(f(3))) = 3;

buildHistogram(imnoisy,configurationTruth);

Etrue = calcLikelihood(miuf,imnoisy,configurationTruth,alpha) + calcSmoothnessPrior(configurationTruth, beta);

%% Q2
% test = [1 0 1 0 0;0 1 1 1 2;1 0 0 3 5;1 0 1 5 1;8 2 3 1 1;];
% V2 = calcSmoothnessPrior(im1, 1);
% V21 = burteForceSol(im1, 1);

figure;
subplot(3,2,1);imshow(imreal); title('Raw','FontName','Arial','FontSize',20);
subplot(3,2,2);imshow(imnoisy); title('Noisy','FontName','Arial','FontSize',20);


% configuration 1
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

% configuration 2
configuration2 = medfilt2(configuration1,[3,3],'symmetric');
imseg2 = segmentation(imnoisy, configuration2, f, miuf);
L2 = calcLikelihood(miuf,imnoisy,configuration2,alpha);
P2 = calcSmoothnessPrior(configuration2, beta);
E2 = L2 + P2;
subplot(3,2,4);imshow(imseg2); title('3x3 Median Filter','FontName','Arial','FontSize',20);

% configuration 3
imsmmoth = imgaussfilt(imnoisy,1,'FilterSize',3,'Padding','replicate');
configuration3(imsmmoth<100) = f(1);
configuration3(imsmmoth>=100 & imsmmoth <= 160) = f(2);
configuration3((imsmmoth>=160)) = f(3);

imseg3 = segmentation(imnoisy, configuration3, f, miuf);
L3 = calcLikelihood(miuf,imnoisy,configuration3,alpha);
P3 = calcSmoothnessPrior(configuration3, beta);
E3 = L3 + P3
subplot(3,2,5);imshow(imseg3); title('Gaussian','FontName','Arial','FontSize',20);

% configuration 4
se = strel('disk',1);
configuration4 = imopen(configuration1,se);
imseg4 = segmentation(imnoisy, configuration4, f, miuf);
L4 = calcLikelihood(miuf,imnoisy,configuration4,alpha);
P4 = calcSmoothnessPrior(configuration4, beta);
E4 = L4 + P4;
subplot(3,2,6);imshow(imseg4); title('Morpho','FontName','Arial','FontSize',20);

figure
subplot(2,2,1);
buildHistogram(imnoisy,configuration1); title('Threshold','FontName','Arial','FontSize',20);
subplot(2,2,2);
buildHistogram(imnoisy,configuration2); title('3x3 Median Filter','FontName','Arial','FontSize',20);
subplot(2,2,3);
buildHistogram(imnoisy,configuration3); title('Gaussian filter','FontName','Arial','FontSize',20);
subplot(2,2,4);
buildHistogram(imnoisy,configuration4); title('Morphological: Close','FontName','Arial','FontSize',20);

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

% try with parallel optimization
iter = 1;
maxiter = 1e2;
oldCost = 1e6;

%configuration = configuration1;% start with cfg1
configuration = randi(3,size(imnoisy));

figure
while iter < maxiter
    configuration = GibbsSampling(imnoisy, configuration, f, miuf, alpha, beta);%ICM
    
    Lnew = calcLikelihood(miuf,imnoisy,configuration,alpha);
    Pnew = calcSmoothnessPrior(configuration, beta);
    newCost = Lnew + Pnew;
    
    if abs(newCost-oldCost) < 1e-3 || newCost > oldCost
        break;
    end
%     configuration = newconfiguration;
    oldCost = newCost;
    costs(iter) = newCost;
    iter = iter + 1;
    imres = segmentation(imnoisy, configuration, f, miuf);
    imshow(imres);pause(0.1);
end
% cost variantion
figure
plot(costs,'r-o','LineWidth',2);grid on;
xlabel('Iteration');
ylabel('Cost')
title('ICM')
set(gca,'FontName','Arial','FontSize',20);

% plot posterior cost
figure
bar([E1 E2 E3 E4 newCost]);
xticklabels({'CF1','CF2','CF3','CF4','ICM'});
xtickangle(45);
set(gca,'FontName','Arial','FontSize',20);
xlim=get(gca,'xlim');
hold on
h1 = plot(xlim,[Etrue Etrue],'r--');
legend([h1],{'Truth'});
title('Posterior');

function configuration = SimulatedAnn
% TODO
end

function configuration = GibbsSampling(im, configuration, f, miuf, alpha, beta)
    mask = checkerboard(1,round(size(im,1)*0.5),round(size(im,2)*0.5));
    mask(size(im,1)+1:end,:) = [];
    mask(:,size(im,2)+1:end) = [];
    mask(mask > 0.5) = 1; mask = mask == 1;
    
    localPotential = calclocalPotentials(im, configuration, f, miuf, alpha, beta);
    prob = exp(-localPotential);
    probNorm = sum(prob,3);
    prob(:,:,1) = prob(:,:,1) ./ probNorm;
    for i = 2:size(localPotential,3) 
        prob(:,:,i) = prob(:,:,i) ./ probNorm;
        prob(:,:,i) = prob(:,:,i) + prob(:,:,i-1);
    end
    
    randProb = rand(size(probNorm));
    
    id = cumsum(prob>randProb,3);
    id = id(:,:,end);
    id = size(localPotential,3) - id;
    id = id + 1;
    
    configuration(mask) = f(id(mask));
    
    localPotential = calclocalPotentials(im, configuration, f, miuf, alpha, beta);
    prob = exp(-localPotential);
    probNorm = sum(prob,3);
    prob(:,:,1) = prob(:,:,1) ./ probNorm;
    for i = 2:size(localPotential,3) 
        prob(:,:,i) = prob(:,:,i) ./ probNorm;
        prob(:,:,i) = prob(:,:,i) + prob(:,:,i-1);
    end
    randProb = rand(size(probNorm));
    id = cumsum(prob>randProb,3);
    id = id(:,:,end);
    id = size(localPotential,3) - id;
    id = id + 1;
    
    configuration(~mask) = f(id(~mask));
end

function imseg = segmentation(im, configuration, f, miuf)
    imseg = im;
    for i = 1:length(f)
        imseg(configuration == f(i)) = miuf(f(i));
    end
end

function V1 = calcLikelihood(miuf, ds, fs, alpha)
    ds = double(ds);
    V1 = alpha*sum((miuf(fs(:))-ds(:)).^2);
end

function V2 = calcSmoothnessPrior(ds, beta)
    ds = double(ds);
    [m,n]=size(ds);
    V2 = 0;
    ds = double(ds);
    % row sweeping
    for i = 1:m-1   
        V2 = V2+sum(abs(ds(i,:)-ds(i+1,:))>1e-3);
    end
    % col sweeping
    for i = 1:n-1   
        V2 = V2+sum(abs(ds(:,i)-ds(:,i+1))>1e-3);
    end
    V2 = V2.*beta;
end

function V2 = burteForceSol(ds,beta)
    [m,n]=size(ds);
    V2 = 0;

    ds = double(ds);
    for i = 1:m
        for j = 1:n
%             if i - 1 > 0
%                 V2 = V2 + double(abs(ds(i,j) - ds(i-1,j))>1e-3);
%             end
            if i + 1 <= m
                V2 = V2 + double(abs(ds(i,j) - ds(i+1,j))>1e-3);
            end
%             if j - 1 > 0 
%                 V2 = V2 + double(abs(ds(i,j) - ds(i,j-1))>1e-3);
%             end
            if j + 1 <= n 
                V2 = V2 + double(abs(ds(i,j) - ds(i,j+1))>1e-3);
            end
        end
    end
    V2 = V2.*beta;
end

function buildHistogram(im1,seg)
    [h,wout] = hist(double(im1(:)),256);
    [h1,wout1] = hist(double(im1(seg==1)),256);  
    [h2,wout2] = hist(double(im1(seg==2)),256);
    [h3,wout3] = hist(double(im1(seg==3)),256);
    
    bar(wout,h);hold on;grid on;

    ccmap = lines(3);
    plot(wout1(h1~=0),h1(h1~=0),'-','LineWidth',2,'Color',ccmap(1,:));
    plot(wout2(h2~=0),h2(h2~=0),'-','LineWidth',2,'Color',ccmap(2,:));
    plot(wout3(h3~=0),h3(h3~=0),'-','LineWidth',2,'Color',ccmap(3,:));
    title('Histogram','FontName','Arial','FontSize',20);
end
