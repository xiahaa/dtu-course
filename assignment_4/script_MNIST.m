clc;close all;clear all;

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

training = 0;

if training == 1
    %% load data
    load('MNIST.mat');
    train_data = single(train_data);
    train_label = single(train_label);

    %% data preprocessing
    % preparation, normalization and minus mean
    [data,mdata,scale] = normalization(train_data);
    % random sample train set and validaton set
    valID = randperm(size(data,1),10000);
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
    opts.train.numEpochs = 500;

    % set learning rate
    opts.train.learningRate = 0.1;
    opts.train.weightDecay  = 1e-4;
    opts.train.momentum = 0.9;
    opts = genBatchIndex(N,opts);
    
    opts.earlyStopping.n = 1;
    opts.earlyStopping.patience = 10;
    opts.earlyStopping.num = 0;
    opts.earlyStopping.bestValLoss = 1e8;
    opts.earlyStopping.nn = {};
    opts.earlyStopping.numEpoches = 0;

    %% define forward neural network
    num_of_hidden_units = [1000 1000];% mxn: m is the hidden units per layer, n - layer 500
    num_inputs = size(trainData,1);
    num_outputs = 10;
    nn = cell(1,length(num_of_hidden_units)+1);
    [nn, opts] = initializeNN(num_of_hidden_units,num_inputs,num_outputs,nn,opts);


    %% training
    figure;
    trainlosses = zeros(1,opts.train.numEpochs);
    vallosses = zeros(1,opts.train.numEpochs);
    % per epoch
    for i = 1:opts.train.numEpochs    
        % gen new batches
        opts = getBatches(N, opts);
        % per batch
        for j = 1:opts.batchNum
            id = opts.batches(j,:);
            x = trainData(:,id);
            t = trainLabel(:,id);
            % forward
            [y,h,z] = forwardPropagation(nn,x,@ReLU);

            newloss = loss(t,y);
            disp(['Training Epoch ',num2str(i),'--|--batch ',num2str(j),':', num2str(newloss)]);

            % backpropagation + Minibatch + SGD + Momentum 0.9
            grad = backpropagation(nn,t,y,h,z,@ReLUDer);        
            % SDG: gradient descent
    %         nn = cellfun(@updateW,nn,grad,opts,'UniformOutput', false);
            nn = updateW(nn,grad,opts);
        end
        % forward
        [y,~,~] = forwardPropagation(nn,trainData,@ReLU);
        % compute loss
        newloss = loss(trainLabel,y);
        disp(['Training Loss Epoch ',num2str(i),': ', num2str(newloss)]);
        
        % check validation error
%         if mod(i,opts.earlyStopping.n) == 0
            [yv,~,~] = forwardPropagation(nn,valData,@ReLU);
            valLoss = loss(valLabel,yv);
            disp(['Vlidation Loss: ', num2str(valLoss)]);
            if valLoss < opts.earlyStopping.bestValLoss
                opts.earlyStopping.num = 0;
                opts.earlyStopping.bestValLoss = valLoss;
                opts.earlyStopping.numEpoches = i;
                opts.earlyStopping.nn = nn;
            else
                opts.earlyStopping.num = opts.earlyStopping.num + 1;
                opts.train.learningRate = opts.train.learningRate * 0.5;
                if opts.train.learningRate < 0.001
                    opts.train.learningRate = 0.001;
                end
            end
            if opts.earlyStopping.num >= opts.earlyStopping.patience
                break;
            end
            vallosses(i) = valLoss;
            
%         end
        
        trainlosses(i) = newloss;
        plot(trainlosses(1:i),'-o','LineWidth',1.5);hold on;grid on;
        plot(vallosses(1:i),'-o','LineWidth',1.5);
        hold off;
        pause(0.1);
    end
    
    %     nn = opts.earlyStopping.nn;
    
    %% combine train data and validation data, re-start a training using found epoches
    trainData = cat(2,trainData,valData);
    trainLabel = cat(2,trainLabel,valLabel);
    N = length(trainData);
    opts = genBatchIndex(N,opts);
    [nn, opts] = initializeNN(num_of_hidden_units,num_inputs,num_outputs,nn,opts);
    opts.train.numEpochs = opts.earlyStopping.numEpoches;
    opts.train.learningRate = 0.1;
    
    %% retraining
    figure;
    trainlosses = zeros(1,opts.train.numEpochs);
    % per epoch
    for i = 1:opts.train.numEpochs    
        % gen new batches
        opts = getBatches(N, opts);
        % per batch
        for j = 1:opts.batchNum
            id = opts.batches(j,:);
            x = trainData(:,id);
            t = trainLabel(:,id);
            % forward
            [y,h,z] = forwardPropagation(nn,x,@ReLU);

            newloss = loss(t,y);
            disp(['Training Epoch ',num2str(i),'--|--batch ',num2str(j),':', num2str(newloss)]);

            % backpropagation + Minibatch + SGD + Momentum 0.9
            grad = backpropagation(nn,t,y,h,z,@ReLUDer);        
            % SDG: gradient descent
    %         nn = cellfun(@updateW,nn,grad,opts,'UniformOutput', false);
            nn = updateW(nn,grad,opts);
        end
        % forward
        [y,~,~] = forwardPropagation(nn,trainData,@ReLU);
        % compute loss
        newloss = loss(trainLabel,y);
        disp(['Training Loss Epoch ',num2str(i),': ', num2str(newloss)]);
        trainlosses(i) = newloss;
        plot(trainlosses(1:i),'-o','LineWidth',1.5);hold on;grid on;
        hold off;
        pause(0.1);
    end  
    
    net.nn = nn;
    net.meandata = mdata;
    net.scale = scale;
    files = dir(fullfile('./mnist_train/', '*.mat'));
    netcnt = length(files) + 1;
    save(strcat('./mnist_train/net',num2str(netcnt),'.mat'),'net');
    
else
    %% load data
    load('MNIST.mat');
    
    files = dir(fullfile('./mnist_train/', '*.mat'));
    load(strcat('./mnist_train/',files(end).name));
    test_data = single(test_data)';
    test_label = single(test_label)';
    test_data = test_data - repmat(net.meandata',1,size(test_data,2));
    test_data = test_data.*net.scale;
%     data = normalization(test_data);
    
    [y,~,~] = forwardPropagation(net.nn,test_data,@ReLU);
    [~,id1]=max(y);
    [~,id2]=max(test_label);
    disp(sum(id1==id2)/length(id2));
    
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
        w = initialization(column, row*column);
        A = reshape(w,row,column);
        nn{i} = A;

        opts.train.v{i} = zeros(size(A));
    end
end

function opts = genBatchIndex(n,opts)
% run once, generate bacth index.
    startIndex = 1;
    i = 1;
    opts.batchIndex = zeros(floor(n/opts.train.batchSize)+1,opts.train.batchSize);
    while startIndex <= n
        endIndex = min((startIndex+opts.train.batchSize-1),n);
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
    m = mean(datan);
    datan = datan - repmat(m,size(data,1),1);
    
    scale = 1;%sqrt(1 ./ mean(diag(datan'*datan)));
    datan = datan .* scale;
end

function W = updateW(W,grad,opts)
    % update velocity
    v = cell(length(W),1);
    for i = 1:length(W)
        v{i} = opts.train.v{i}.*opts.train.momentum - grad{i}.*opts.train.learningRate;%.*(1-opts.train.momentum);
        W{i} = W{i} + v{i};
        opts.train.v{i} = v{i};
    end
end

function w = initialization(n, m)
% n is assumed to be the number of data. m is assumed to be the number of
% weights for the forward neural network

    % draw from a gaussian with sqrt(2/n) as the standard deviation
    w = randn([m,1]).*sqrt(2/n);
end

function h = ReLU(z)
    h = z;
    h(h<0) = 0;
end

function y = softMax(yhat)
    ey = exp(yhat);
    y = ey ./ (sum(ey)+1e-8);
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
    
    h{1} = [x];
    for i = 1:length(W)-1
        % layer to layer
        z{i} = W{i}*[h{i};ones(1,size(h{i},2))];
        h{i+1} = f(z{i});
    end
    % output layer
    z{i+1} = W{i+1}*[h{i+1};ones(1,size(h{i+1},2))];
    y = softMax(z{i+1});
end

function grad = backpropagation(W,t,y,h,z,fgrad)
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
    grad{n} = A2*A1;
    grad{n} = scale .* grad{n};
    for i = length(W)-1:-1:1
        d1 = W{i+1}'*delta{i+1};
        delta{i} = fgrad(z{i}).*d1(1:end-1,:);
        grad{i} = delta{i}*[h{i}' ones(size(h{i},2),1)];
        grad{i} = scale .* grad{i};
    end
end

function flag = isEarlyStopping(trnErrs, valErrs)
% TODO: implementation the early stopping for the NN.
end