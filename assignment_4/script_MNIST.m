clc;close all;clear all;

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

% load
load('MNIST.mat');
train_data = single(train_data);
train_label = single(train_label);

% preparation, normalization and minus mean
data = normalization(train_data);

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

% initialize parameters
opts.train.batchSize = 100;
opts.train.numEpochs = 20;

% set learning rate
opts.train.learningRate = 1e-3;
opts.train.weightDecay  = 1e-4;
opts.train.momentum = 0.9;
opts.train.v = 0;
opts = genBatchIndex(N,opts);

%% define forward neural network
num_of_hidden_units = [1000 300];% mxn: m is the hidden units per layer, n - layer
num_inputs = size(trainData,1);
num_outputs = 10;

htb = zeros(1,size(num_of_hidden_units,1)+1);
ztb = zeros(1,size(num_of_hidden_units,1)+1);

for i = 1:size(num_of_hidden_units,1)+1
    if i == 1 
        % first layer: inputs+bias
        row = num_of_hidden_units(i);
        column = (num_inputs + 1);
    elseif i < size(num_of_hidden_units,1)+1
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
    
    htb(i) = column-1;
    ztb(i) = row;
end

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
        % compute loss
        newloss = loss(t,y);
        % SGD: backpropagation
        hc = cell2mat(h);
        zc = cell2mat(z);
        grad = cell(opts.train.batchSize,1)';
        for k = 1:size(y,2)
            hk = mat2cell(hc(:,k),htb);
            zk = mat2cell(zc(:,k),ztb);
            grad{k} = backpropagation(nn,t(:,k),y(:,k),hk,zk,@ReLUDer);
        end
        
        % SDG: gradient descent
        nn = cellfun(@updateW,nn,grad,'UniformOutput', false);

    %     cyclic_id = cyclic_id + 1;if cyclic_id > datalength; cyclic_id = 1;end
        cyclic_id = randperm(datalength,1);

        iter = iter + 1;
        losses(iter) = newloss;
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

function datan = normalization(data)
% normalization for MNIST
    datan = data ./ vecnorm(data,2,2);
    m = mean(data);
    datan = datan - repmat(m,size(data,1),1);
end

function W = updateW(W,grad,opts)
    % update velocity
    v = opts.train.v.*opts.train.momentum - grad.*opts.train.learningRate;
    W = W + v;
end

function w = initialization(n, m)
% n is assumed to be the number of data. m is assumed to be the number of
% weights for the forward neural network

    % draw from a gaussian with sqrt(2/n) as the standard deviation
    w = rand([m,1]).*sqrt(2/n);
end

function h = ReLU(z)
    h = z;
    h(h<0) = 0;
end

function y = softMax(yhat)
    ey = exp(yhat);
    y = ey ./ sum(ey);
end

function hdz = ReLUDer(z)
    hdz = ones(numel(z),1);
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
    n = length(W);
    delta{n} =  y - t;
    A1 = repmat([h{n}' 1], numel(delta{n}), 1);
    A2 = repmat(delta{n}, 1, size(h{n},1)+1);
    grad{n} = A1.*A2;
    for i = length(W)-1:1
        d1 = W{i+1}'*delta{i+1};
        delta{i} = fgrad(z{i}).*d1(1:end-1);
        A1 = repmat([h{i}' 1], numel(delta{i}), 1);
        A2 = repmat(delta{i}, 1, size(h{i},1)+1);
        grad{i} = A1.*A2;
    end
end

function flag = isEarlyStopping(trnErrs, valErrs)
% TODO: implementation the early stopping for the NN.
end