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
trainCNN = 1;

id = [16 18 20 21];
nnnum = 1:20;% test nearest neightbors

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
elseif trainCNN == 0
    %% only use the extracted features for classification
    doTraining = 0;
    
    if doTraining == 1
        for i = 1:length(id)
            % load features
            f1 = load(['features' num2str(id(i)) '.mat']);
            % get structure fields
            fields = fieldnames(f1);
            % rename to a uniform name
            eval(['feature{',num2str(i),'}=f1.',fields{1},';']);
        end
        % per features and per nnnum
        errors = zeros(length(feature),length(nnnum));

        % now, we have several feature vectors, we use knn as the classfier.
        for i = 1:length(feature)
            fea = feature{i};
            parfor j = 1:length(nnnum)
                errors(i,j) = trainKNN(fea, nnnum(j), labels);
            end
        end
        save('errors.mat','errors');
    else
        load('errors.mat');
    end
    
    % plot
    figure;
    ccmap = lines(size(errors,1));
    for i = 1:size(errors,1)
        subplot(2,2,i);
        plot(nnnum,errors(i,:),'-d','LineWidth',1.5,'Color',ccmap(i,:));
        xticks(nnnum);
        xtickangle(45);
        title(['Feature Layer: ',num2str(id(i))]);
        xlabel('k neighbors');
        ylabel('accuracy');
        set(gca,'FontName','Arial','FontSize',15);
        grid on;
    end
    
    % best per feature
    [bestScores,maxid] = max(errors,[],2);
    figure;
    bar(bestScores); grid on;
    tickslabel = cell(1,length(id));
    for i = 1:length(id)
        tickslabel{i} = strcat('Layer ',num2str(id(i)),'-',num2str(nnnum(maxid(i))));
    end
    xticklabels(tickslabel);
    xtickangle(45);
    title(['Feature Comparison']);
    xlabel('Feature layers');
    ylabel('Best accuracy');
    set(gca,'FontName','Arial','FontSize',15);
else
    training = 0;
    
    if training == 1
        %% here, I tried to train only the last several layers of the pretrained mode to fit new data
        load('features15.mat');

        trainNum = round(size(x,4) * 0.5);

        imdb.images.data = x;
        imdb.images.mean = zeros(size(x,1),size(x,2),size(x,3));% do not do normalization
        imdb.images.label = single(labels');
        imdb.images.set = [ones(1,trainNum) 3*ones(1,round(size(x,4)*0.6)-trainNum)];% 50% as the training set, 10% as the validation set
        imdb.meta.sets = {'train','val','test'};

        % initialize parameters
        opts.train.batchSize = 100;
        opts.train.numEpochs = 20;

        opts.train.continue = false;
        % use gpu
        opts.train.gpus = [];

        % set learning rate
        opts.train.learningRate = 1e-3;

        opts.train.expDir = 'epoch_data_1';
        opts.train.numSubBatches = 1;

        opts.train.weightDecay  = 1e-4;

        % define a new subnet
        f=1 / 100 ;
        net.layers = {};
        %Fully Connected Layer
        net.layers{end+1} = struct('name', 'fc6', ...
                                   'type', 'conv', ...
                                   'weights', {{f*randn(6,6,256,4096, 'single'), zeros(4096, 1, 'single')}}, ...
                                   'stride', [1, 1], ...
                                   'size', [6,6,256,4096], ...
                                   'pad', [0, 0, 0, 0]);
        %ReLU Layer
        net.layers{end+1} = struct('name','relu6','type', 'relu') ;
        %Fully Connected Layer
        net.layers{end+1} = struct('name', 'fc7', ...
                                   'type', 'conv', ...
                                   'weights', {{f*randn(1,1,4096,4096, 'single'), zeros(4096, 1, 'single')}}, ...
                                   'stride', [1, 1], ...
                                   'size', [1,1,4096,4096], ...
                                   'pad', [0, 0, 0, 0]) ;
        %ReLU Layer
        net.layers{end+1} = struct('name','relu7','type', 'relu') ;
        %Fully Connected Layer
        net.layers{end+1} = struct('name', 'fc8', ...
                                   'type', 'conv', ...
                                   'weights', {{f*randn(1,1,4096,2, 'single'), zeros(2, 1, 'single')}}, ...
                                   'stride', [1, 1], ...
                                   'size', [1,1,4096,2], ...
                                   'pad', [0, 0, 0, 0]) ;
        %Output Layer: softmaxloss
        net.layers{end+1} = struct('name','objective','type', 'softmaxloss') ;
        net = vl_simplenn_tidy(net) ;

        % Convert to a GPU array if needed
        if numel(opts.train.gpus) > 0
            imdb.images.data = gpuArray(imdb.images.data) ;
        end

        % Call training function in MatConvNet
        [net,info] = cnn_train(net, imdb, @getBatch, opts.train,'val',find(imdb.images.set==3)) ;

        % Move the CNN back to the CPU if it was trained on the GPU
        if numel(opts.train.gpus) > 0
            net = vl_simplenn_move(net, 'cpu') ;
        end

        save('data/retraincnn_new.mat', '-struct', 'net') ;
    else
        if 0
            net1 = load('data/retraincnn_new.mat');
            net2 = load(strcat(data_path,'imagenet-vgg-f.mat'));

            net.layers = {};
            for i = 1:15
                net.layers{end+1} = net2.layers{i};
            end
            for i = 1:5
                net.layers{end+1} = net1.layers{i};
            end
            net.layers{end+1} = struct('name','prob','type', 'softmax') ;
            net = vl_simplenn_tidy(net);

            % preprocess input images
            im_ = cnnPreprocess(histo_images, net2);
            % Run the CNN.
            res = vl_simplenn(net, im_) ;

             % Show the classification result, just debug
            scores = squeeze(gather(res(end).x)) ;
            [bestScore, best] = max(scores) ;

            % save result
            save('cnnres.mat','best');
        end
        
        %% finanlly, compare CNN with KNN
        load('errors.mat');
        % best per feature
        [bestScores,maxid] = max(errors,[],2);
        
        % CNN
        load('cnnres.mat');
        correctAns = sum(uint8(best') == labels)/length(labels);
        
        figure;
        bar([bestScores;correctAns]); grid on;
        tickslabel = cell(1,length(id));
        for i = 1:length(id)
            tickslabel{i} = strcat('Layer ',num2str(id(i)),'-',num2str(nnnum(maxid(i))));
        end
        tickslabel{end+1} = 'CNN retrained';        
        xticklabels(tickslabel);
        xtickangle(45);
        title(['Feature Comparison']);
        xlabel('Feature layers');
        ylabel('Best accuracy');
        set(gca,'FontName','Arial','FontSize',15);
    end
    
end

% --------------------------------------------------------------------
function [images,labels] = getBatch(imdb, batch)
    % --------------------------------------------------------------------
    images = imdb.images.data(:,:,:,batch) ;
    labels = imdb.images.label(:,batch) ;
    
%     if opts.train.gpus > 0
%   		images = gpuArray(images);
%     end
%     inputs = {'input', images, 'label', labels} ;
end

function varargout = KNNClassifier(train, query, nn, labels)
    % use Euclidean distance
    eucDist = pdist2(train',query');
    % find nearest neighbors
    [eucDistSorted,orders] = sort(eucDist,'ascend');

    % find most support label as the label of the test sample
    testLabel = mode(labels(orders(1:nn)));
    
    varargout{1} = testLabel;
    varargout{2} = eucDistSorted(1:nn);
end

%% todo: write a knn classifier and use leave-one out validation to analyze which combination is the best.
% for every feature, different kind of knn classifier will be compared.
% compare the best knn for different features.
function varargout = trainKNN(feature, nn, labels)
    % sample num
    n = size(feature,2);
    % generate random permuted seeds for validation
    seeds = randperm(n,n);
    tt = 1:n;
    cnum = 0;
    for i = 1:length(seeds)
        % get train id and test id
        testid = seeds(i);
        trainid = tt(tt~=testid);
        % data
        testData = feature(:,testid);
        trainData = feature(:,trainid);
        
        % KNN classifier
        testLabel = KNNClassifier(trainData, testData, nn, labels(trainid));
        
        % check if test label is correct
        if testLabel == labels(testid)
            cnum = cnum + 1;
        end
    end
    % leave-one-out error
    varargout{1} = cnum ./ n;
end

function im = cnnPreprocess(batch_images, net)
    % Preprocess images
    im = single(batch_images);
    im = imresize(im, net.meta.normalization.imageSize(1:2));
	im = bsxfun(@minus,im,net.meta.normalization.averageImage);
end

