clc;close all;clear all;

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

%% load data
load('MNIST.mat');
train_data = single(train_data);
train_label = single(train_label);

test_data = single(test_data);
test_label = single(test_label);

% load('imdb.mat');

%%
im = reshape(train_data',28,28,1,[]);
% im = permute(im,[2,1,3,4]);

useBuiltinFunc = true;

%%%%%%%%%%%%%%%%%%%%%%% injecting noise, scaling, rotation and translation %%%%%%%%%%%%%%%%%%%%%%%%%%
if useBuiltinFunc == false
    [~,labels] = max(train_label');
    imns = [];
    lbns = [];
    for i = 1:10
        id = find(labels == i);
        % 
        idn = round(length(id) * 0.25);
        ids = id(randperm(length(id),idn));

        imn = zeros(28,28,1,length(ids)*4);

        for j = 1:1:length(ids)
            imn(:,:,:,(j-1)*4+1) = randomNoise(im(:,:,:,ids(j)),0.3);
            imn(:,:,:,(j-1)*4+2) = rescaling(im(:,:,:,ids(j)));
            imn(:,:,:,(j-1)*4+3) = randomRot(im(:,:,:,ids(j)),10);
            imn(:,:,:,(j-1)*4+4) = randomTrans(im(:,:,:,ids(j)),5);
        end
        imns = cat(4,imns,imn);
        lb = single(zeros(10,size(imn,4)));
        lb(i,:) = 1;
        lbns = cat(2,lbns,lb);
    end
    imns = reshape(imns,28*28,[])';
    save('aug3.mat','imns','lbns');
else
    images = reshape(train_data',28,28,1,[]);
    [~,labels] = max(train_label');
    
    augmenter = imageDataAugmenter( ...
                                    'RandRotation',[-20 20], ...
                                    'RandXScale',[0.6 1.1], ...
                                    'RandYScale',[0.6 1.1], ...
                                    'RandXTranslation',[-3 3], ...
                                    'RandYTranslation',[-3 3]);
                                
    imageSize = [28 28 1];
    augimds = augmentedImageDatastore(imageSize,images,labels(:),'DataAugmentation',augmenter);
    
    num = floor(length(labels)/128);
    
    aug_labels = zeros(num*128,10);
    for i = 1:num
        minibatch = augimds.read();
        montage(minibatch.input);
        pause(0.1);
        images = cat(4,images,minibatch.input{:});
        for j = 1:length(minibatch.response)
            aug_labels((i-1)*128+j, minibatch.response{j}) = 1; 
        end
    end
    
    aug_data = reshape(images,784,[])';
    aug_label = cat(1,train_label,aug_labels);
%     minibatch = preview(augimds);
%     montage(minibatch.input);
%     
    save('aug_data4.mat','aug_data','aug_label');
end

function imt = randomRot(im,maxdeg)
% max 10, ok
    ang = (rand(1)*2 - 1) * pi / 180.0;%(rand(1)*2 - 1)*
    tform = affine2d([cos(ang) sin(ang) 0;-sin(ang) cos(ang) 0;0 0 1]');
    R = imref2d(size(im));
    imt = imwarp(im,tform,'InterpolationMethod','bilinear','OutputView',R);
end

function imt = randomTrans(im,maxT)
% max 5, ok
    t = (rand(2,1).*2-1).*maxT;%
    tform = affine2d([1 0 t(1);0 1 t(2);0 0 1]');
    R = imref2d(size(im));
    imt = imwarp(im,tform,'InterpolationMethod','bilinear','OutputView',R);
end

function imt = randomNoise(im,sigma)
% 0.2 = sigma
    imt = im + randn([size(im,1),size(im,2)]).*sigma;
end

function imt = rescaling(im)
% function like smoothing
%     ims = imresize(im,0.8,'bilinear');
%     imt = imresize(ims,size(im),'bilinear');
    t = [1;1]+rand([2,1]).*0.6-0.3;
    tform = affine2d([t(1) 0 0;0 t(2) 0;0 0 1]');
    R = imref2d(size(im));
    imt = imwarp(im,tform,'InterpolationMethod','bilinear','OutputView',R);
end
