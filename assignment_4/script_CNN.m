clc;close all;clear all;

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

% load data
data_path = '../data/Ex_6_data/';

% 4d data: rows x columns x channels x number of frames
load(strcat(data_path,'histo_images.mat'));

% visualization
num = size(histo_images,4);

UseCNN = 0;

id = [16 18 20 21];

if UseCNN == 1
    %% forward section using pretrained model
    net = load(strcat(data_path,'imagenet-vgg-f.mat'));
    net = vl_simplenn_tidy(net);
    
    % preprocess input images
    im_ = cnnPreprocess(histo_images, net);
    % Run the CNN.
    res = vl_simplenn(net, im_) ;

    % Show the classification result, just debug
    scores = squeeze(gather(res(end).x)) ;
    [bestScore, best] = max(scores) ;

    % for i = 1:num
    %     im = histo_images(:,:,:,i);
    %     im_ = single(im);
    %     % preprocessing
    %     im_ = imresize(im_,net.meta.normalization.imageSize(1:2));
    %     im_ = im_ - net.meta.normalization.averageImage;
    %     
    %     % Run the CNN.
    %     res = vl_simplenn(net, im_) ;
    %     
    %     % Show the classification result.
    %     scores = squeeze(gather(res(end).x)) ;
    %     [bestScore, best] = max(scores) ;
    %     figure(1) ; clf ; imagesc(im) ; 
    %     title(sprintf('%s (%d), score %.3f',...
    %     net.meta.classes.description{best}, best, bestScore)) ;    
    %     pause(0.1);
    % end
    
    % layers id
    for i = 1:length(id)
        eval(['features' num2str(id(i)) '= squeeze(gather(res(id(i)+1).x));']);
        save(['features' num2str(id(i)) '.mat'],['features' num2str(id(i))]);
    end
else
    %% only use the extracted features for classification
    nnnum = [1 2 3 4 5 6 7 8 9 10];% test nearest neightbors
    
    for i = 1:length(id)
        load(['features' num2str(id(i)) '.mat']);
        % now, we have several feature vectors, we use knn as the classfier.
        
        
    end
    
    
    
end




function im = cnnPreprocess(batch_images, net)
    % Preprocess images
    im = single(batch_images);
    im = imresize(im, net.meta.normalization.imageSize(1:2));
	im = bsxfun(@minus,im,net.meta.normalization.averageImage);
end

