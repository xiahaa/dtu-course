clc;close all;clear all;

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

addpath ../utils;

data_dir = '../data/EX_3_data/';

imnoisy = imread(strcat(data_dir,'noisy_circles.png'));
imreal = imread(strcat(data_dir,'noise_free_circles.png'));

f = [1;2;3];
miuf = [70;130;190];

%% Q1
configurationTruth = imreal;
configurationTruth(configurationTruth == miuf(f(1))) = 1;
configurationTruth(configurationTruth == miuf(f(2))) = 2;
configurationTruth(configurationTruth == miuf(f(3))) = 3;

buildHistogram(imnoisy,configurationTruth);

%% Q2
% test = [1 0 1 0 0;0 1 1 1 2;1 0 0 3 5;1 0 1 5 1;8 2 3 1 1;];
% V2 = calcSmoothnessPrior(im1, 1);
% V21 = burteForceSol(im1, 1);

alpha = 0.0005;

figure;
subplot(3,2,1);imshow(imreal); title('Raw','FontName','Arial','FontSize',20);
subplot(3,2,2);imshow(imnoisy); title('Noisy','FontName','Arial','FontSize',20);


%% configuration 1
threshold1 = 100;
threshold2 = 160;
configuration1 = imnoisy;
configuration1(imnoisy<100) = f(1);
configuration1(imnoisy>=100 & imnoisy <= 160) = f(2);
configuration1((imnoisy>=160)) = f(3);

imseg1 = segmentation(imnoisy, configuration1, f, miuf);
E1 = calcLikelihood(miuf,imnoisy,configuration1,alpha) + calcSmoothnessPrior(configuration1, 1);
subplot(3,2,3);imshow(imseg1); title('Threshold','FontName','Arial','FontSize',20);

%% configuration 2
configuration2 = medfilt2(configuration1,[3,3],'symmetric');
imseg2 = segmentation(imnoisy, configuration2, f, miuf);
E2 = calcLikelihood(miuf,imnoisy,configuration2,alpha) + calcSmoothnessPrior(configuration2, 1);
subplot(3,2,4);imshow(imseg2); title('3x3 Median Filter','FontName','Arial','FontSize',20);

%% configuration 3
configuration2 = medfilt2(configuration1,[3,3],'symmetric');
imseg2 = segmentation(imnoisy, configuration2, f, miuf);
E2 = calcLikelihood(miuf,imnoisy,configuration2,alpha) + calcSmoothnessPrior(configuration2, 1);
subplot(3,2,4);imshow(imseg2); title('3x3 Median Filter','FontName','Arial','FontSize',20);

%% configuration 4
configuration2 = medfilt2(configuration1,[3,3],'symmetric');
imseg2 = segmentation(imnoisy, configuration2, f, miuf);
E2 = calcLikelihood(miuf,imnoisy,configuration2,alpha) + calcSmoothnessPrior(configuration2, 1);
subplot(3,2,4);imshow(imseg2); title('3x3 Median Filter','FontName','Arial','FontSize',20);

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
