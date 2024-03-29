clc;close all;clear all;

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

training = 0;

if training == 1
    %% load data
    load('MNIST.mat');
    useAugData = true;
    
    if useAugData == false
        train_data = single(train_data);
        train_label = single(train_label);
    else
        load('aug_data8.mat');
        train_data = single(aug_data);
        train_label = single(aug_label);
    end
    
    Nt = size(train_data,1);

    %% data preprocessing
    % preparation, normalization and minus mean
    [data,mdata,scale] = normalization(train_data);
    % random sample train set and validaton set
    valID = randperm(size(data,1),round(Nt*0.1));%% 10% as the validation set
    % loginal id
    id = zeros(1,size(data,1));
    id(valID) = 1;
    id = id == 1;
    % separate data to training and validation
    valData = data(id,:)';
    valLabel = train_label(id,:)';
    trainData = data(~id,:)';
    trainLabel = train_label(~id,:)';
    
    N = length(trainData);

    %% training parameters
    % initialize parameters
    opts.train.batchSize = 100;
    opts.train.numEpochs = 100;

    lrinit = 0.01;
    
    % set learning rate
    opts.train.learningRate = lrinit;
    opts.train.weightDecay  = 0;%1e-4;
    opts.train.momentum = 0.9;
    opts = genBatchIndex(N,opts);
    
    opts.earlyStopping.n = 1;
    opts.earlyStopping.patience = 20;
    opts.earlyStopping.num = 0;
    opts.earlyStopping.bestValLoss = 1e8;
    opts.earlyStopping.nn = {};
    opts.earlyStopping.numEpoches = 0;
    
    %% Adam parameters, if not Adam, then useless
    opts.train.Adam.eps = 0.001;
    opts.train.Adam.rho1 = 0.9;
    opts.train.Adam.rho2 = 0.999;
    opts.train.Adam.delta = 1e-8;
    opts.train.Adam.t = 0;

    %% define forward neural network
    num_of_hidden_units = [1000];% mxn: m is the hidden units per layer, n - layer 1000 300 best
    num_inputs = size(trainData,1);
    num_outputs = 10;
    nn = cell(1,length(num_of_hidden_units)+1);
    [nn, opts] = initializeNN(num_of_hidden_units,num_inputs,num_outputs,nn,opts);
    
    %% just for showing test error
    test_data = single(test_data);
    test_label = single(test_label)';
    test_data = test_data';
    test_data = test_data - repmat(mdata',1,size(test_data,2));
    test_data = test_data.*scale;
    
    %% training
    figure;
    trainlosses = zeros(1,opts.train.numEpochs);
    vallosses = zeros(1,opts.train.numEpochs);
    testlosses = zeros(1,opts.train.numEpochs);
    % per epoch
    for i = 1:opts.train.numEpochs    
        % gen new batches
        opts = getBatches(N, opts);
        % per batch
        newloss = 0;
        for j = 1:opts.batchNum
            id = opts.batches(j,:);
            x = trainData(:,id);
            t = trainLabel(:,id);
            % forward
            [y,h,z] = forwardPropagation(nn,x,@ReLU);

            newlossminibatch = loss(t,y);
            disp(['Training Epoch ',num2str(i),'--|--batch ',num2str(j),': ', num2str(newlossminibatch/length(t))]);

            % backpropagation + Minibatch + SGD + Momentum 0.9
            grad = backpropagation(nn,t,y,h,z,@ReLUDer,opts.train.weightDecay);        
            % SDG: gradient descent
            [nn,opts] = updateW(nn,grad,opts);
            % use Adam
%             opts.train.Adam.t = opts.train.Adam.t + 1;
%             [nn,opts] = updateWAdam(nn,grad,opts);
            newloss = newloss + newlossminibatch;
        end
        
        disp(['Training Loss Epoch ',num2str(i),': ', num2str(newloss/N)]);
        
        % check validation error
        [yv,~,~] = forwardPropagation(nn,valData,@ReLU);
        valLoss = loss(valLabel,yv);
        disp(['Vlidation Loss: ', num2str(valLoss/length(valLabel))]);
        if valLoss < opts.earlyStopping.bestValLoss
            opts.earlyStopping.num = 0;
            opts.earlyStopping.bestValLoss = valLoss;
            opts.earlyStopping.numEpoches = i;
            opts.earlyStopping.nn = nn;
        else
            opts.earlyStopping.num = opts.earlyStopping.num + 1;
        end
        if opts.earlyStopping.num >= opts.earlyStopping.patience
            break;
        end
        vallosses(i) = valLoss/length(valLabel);
        
        % check test accuracy
        [yt,~,~] = forwardPropagation(nn,test_data,@ReLU);
        testLoss = loss(test_label,yt);
        [~,id1]=max(yt);
        [~,id2]=max(test_label);
        disp(['test accuracy: ', num2str(sum(id1==id2)/length(id2))]);
        testlosses(i) = testLoss/length(yt);
        
        % plot
        trainlosses(i) = newloss/N;
        plot(trainlosses(1:i),'-o','LineWidth',1.5);hold on;grid on;
        plot(vallosses(1:i),'-o','LineWidth',1.5);
        plot(testlosses(1:i),'-o','LineWidth',1.5);
        hold off;
        pause(0.1);
    end
    
    bestnn = opts.earlyStopping.nn;
    
    %% combine train data and validation data, re-start a training using found epoches
    trainData = cat(2,trainData,valData);
    trainLabel = cat(2,trainLabel,valLabel);
    N = length(trainData);
    opts = genBatchIndex(N,opts);
    [nn, opts] = initializeNN(num_of_hidden_units,num_inputs,num_outputs,nn,opts);
    opts.train.numEpochs = opts.earlyStopping.numEpoches;
    opts.train.learningRate = lrinit;
    opts.train.Adam.t = 0;
    
    %% retraining
    figure;
    trainlosses = zeros(1,opts.train.numEpochs);
    % per epoch
    for i = 1:opts.train.numEpochs
        % gen new batches
        opts = getBatches(N, opts);
        % per batch
        newloss = 0;
        for j = 1:opts.batchNum
            id = opts.batches(j,:);
            x = trainData(:,id);
            t = trainLabel(:,id);
            % forward
            [y,h,z] = forwardPropagation(nn,x,@ReLU);

            newlossminibatch = loss(t,y);
            
            disp(['Training Epoch ',num2str(i),'--|--batch ',num2str(j),': ', num2str(newlossminibatch/length(t))]);

            % backpropagation + Minibatch + SGD + Momentum 0.9
            grad = backpropagation(nn,t,y,h,z,@ReLUDer,opts.train.weightDecay);        
            % SDG: gradient descent
            [nn,opts] = updateW(nn,grad,opts);
%             opts.train.Adam.t = opts.train.Adam.t + 1;
%             [nn,opts] = updateWAdam(nn,grad,opts);
            
            newloss = newloss + newlossminibatch;
        end
        
        % check test accuracy
        [yt,~,~] = forwardPropagation(nn,test_data,@ReLU);
        testLoss = loss(test_label,yt);
        [~,id1]=max(yt);
        [~,id2]=max(test_label);
        disp(['test accuracy: ', num2str(sum(id1==id2)/length(id2))]);
        testlosses(i) = testLoss/length(yt);

        % plot
        disp(['Training Loss Epoch ',num2str(i),': ', num2str(newloss/N)]);
        trainlosses(i) = newloss/N;
        plot(trainlosses(1:i),'-o','LineWidth',1.5);hold on;grid on;
        plot(testlosses(1:i),'-o','LineWidth',1.5);
        hold off;
        pause(0.1);
    end  
    
    % save net work
    net.nn = nn;
    net.meandata = mdata;
    net.scale = scale;
    files = dir(fullfile('./mnist_train/', '*.mat'));
    netcnt = length(files) + 1;
    save(strcat('./mnist_train/net',num2str(netcnt),'.mat'),'net');
    
else
    %% load data
    load('MNIST.mat');
    
    %% load net
    load(strcat('./mnist_train/','net36.mat'));
    
    %% preprocessing
    test_data = single(test_data);
    test_label = single(test_label)';
    
    test_data = test_data';
    test_data = test_data - repmat(net.meandata',1,size(test_data,2));
    test_data = test_data.*net.scale;

    % apply net on test data
    [y,~,~] = forwardPropagation(net.nn,test_data,@ReLU);
    [~,id1]=max(y);
    [~,id2]=max(test_label);
    
    % compute accuracy
    disp(sum(id1==id2)/length(id2));
    
    % show false classified samples
    id = id1 ~= id2;
    falseData = test_data(:,id);
    falseData = falseData + repmat(net.meandata',1,size(falseData,2));
    falseImages = reshape(falseData,28,28,[]);
    falseImages = permute(falseImages,[2,1,4,3]);
    montage(falseImages);
end

function [nn, opts] = initializeNN(num_of_hidden_units,num_inputs,num_outputs,nn,opts)
    for i = 1:length(num_of_hidden_units)+1
        if i == 1 
            % first layer: inputs+bias
            row = num_of_hidden_units(i);
            column = (num_inputs + 1);
        elseif i < length(num_of_hidden_units)+1
            % intermediate layer: previous hidden units+bias
            row = num_of_hidden_units(i);
            column = (num_of_hidden_units(i-1) + 1);
        else
            % output layer:
            row = num_outputs;
            column = (num_of_hidden_units(i-1) + 1);
        end
        w = initialization(row, column, row*column);
        A = reshape(w,row,column);
        nn{i} = A;

        opts.train.s{i} = zeros(size(A));
        opts.train.r{i} = zeros(size(A));
    end
end

function opts = genBatchIndex(n,opts)
% run once, generate bacth index.
    startIndex = 1;
    i = 1;
    opts.batchIndex = zeros(floor(n/opts.train.batchSize)+1,opts.train.batchSize);
    while startIndex <= n
        endIndex = min((startIndex+opts.train.batchSize-1),n);
        if (endIndex - startIndex + 1) ~= opts.train.batchSize
            break;
        end
        opts.batchIndex(i,:) = startIndex:endIndex;
        i = i + 1;
        startIndex = endIndex + 1;
    end
    opts.batchIndex(i:end,:) = [];
    opts.batches = zeros(size(opts.batchIndex));
    opts.batchNum = i-1;
end

function opts = getBatches(n, opts)
% run per epoch, draw new bacthes.
    rn = randperm(n,n)';
    for i = 1:opts.batchNum
        opts.batches(i,:) = rn(opts.batchIndex(i,:));
    end
end

function res = fvecnorm(x,p,dir)
    if dir == 1
        x = x';
    end
    res = zeros(size(x,1),1);
    for i = 1:size(x,1)
        res(i) = norm(x(1,:),p);
    end
end

function [datan,m,scale] = normalization(data)
% normalization for MNIST
%     datan = data ./ fvecnorm(data,2,2);
    datan = data;
%     datan = datan.*2-1;
    m = mean(datan);
    datan = datan - repmat(m,size(data,1),1);
    scale = 1;%sqrt(1 ./ mean(diag(datan'*datan)));
    datan = datan .* scale;
end

function [W,opts] = updateW(W,grad,opts)
% SGD with momentum
    % update velocity
    v = cell(length(W),1);
    for i = 1:length(W)
        v{i} = opts.train.s{i}.*opts.train.momentum - grad{i}.*opts.train.learningRate;%.*(1-opts.train.momentum);
        W{i} = W{i} + v{i};
        opts.train.s{i} = v{i};
    end
end

function [W,opts] = updateWAdam(W,grad,opts)
% Adam optimization
    c1 = 1 / (1-opts.train.Adam.rho1^opts.train.Adam.t);
    c2 = 1 / (1-opts.train.Adam.rho2^opts.train.Adam.t);
    for i = 1:length(W)
        opts.train.s{i} = opts.train.s{i}.*opts.train.Adam.rho1 + (1-opts.train.Adam.rho1).*grad{i};
        opts.train.r{i} = opts.train.r{i}.*opts.train.Adam.rho2 + (1-opts.train.Adam.rho2).*(grad{i}.*grad{i});
        sb = opts.train.s{i}.*c1;
        rb = opts.train.r{i}.*c2;
        deltaW = -opts.train.Adam.eps.*(sb./(sqrt(rb)+opts.train.Adam.delta));
        W{i} = W{i} + deltaW;
    end
end

function w = initialization(m, n, mn)
% n is assumed to be the number of data. m is assumed to be the number of
% weights for the forward neural network

%     w = randn([mn,1]).*sqrt(2/mn);% draw from a gaussian with sqrt(2/n) as the standard deviation
    w = (rand([mn,1]).*2 - 1).*sqrt(6/(m+n));% draw from uniform distribution
end

function h = ReLU(z)
    h = z;
    h(h<0) = 0;
end

function y = softMax(yhat)
    ey = exp(double(yhat));
    y = ey ./ (sum(ey)+1e-20);
    y = y + 1e-20;
end

function hdz = ReLUDer(z)
    hdz = ones(size(z));
    hdz(z<0) = 0;
end

function L = loss(t,y)
% summation of summation along each class. 
    L = sum(-sum(t.*log(y)));
end

function [y,h,z] = forwardPropagation(W,x,f)
% run forward propagation for the forward nn. W is assumed to be a cell
% with each element being the weight matrix of layer i. f is assumed to be
% the activation function handler.
    z = cell(length(W),1);
    h = cell(length(W),1);
    
    sigma = 0;
    h{1} = [x];
    h{1} = h{1} + randn([size(h{1})]).*sigma;
    for i = 1:length(W)-1
        % layer to layer
        z{i} = W{i}*[h{i};ones(1,size(h{i},2))];
        h{i+1} = f(z{i});
        h{i+1} = h{i+1} + randn([size(h{i+1})]).*sigma;
    end
    % output layer
    z{i+1} = W{i+1}*[h{i+1};ones(1,size(h{i+1},2))];
    y = softMax(z{i+1});
end

function grad = backpropagation(W,t,y,h,z,fgrad,weightDecay)
% run backpropagation for gradient computation. W is assumed to be a cell.
    % first, take care of output
    grad = cell(length(W),1)';
    delta = cell(length(W),1)';
    
    miniBatchSize = size(y,2);
    scale = 1./miniBatchSize;
    
    n = length(W);
    delta{n} =  y - t;
    A1 = [h{n}' ones(size(y,2),1)];
    A2 = delta{n};
    grad{n} = A2*A1 + weightDecay.*W{n};
    grad{n} = scale .* grad{n};
    for i = length(W)-1:-1:1
        d1 = W{i+1}'*delta{i+1};
        delta{i} = fgrad(z{i}).*d1(1:end-1,:);
        grad{i} = delta{i}*[h{i}' ones(size(h{i},2),1)] + weightDecay.*W{i};
        grad{i} = scale .* grad{i};
    end
end
