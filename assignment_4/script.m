clc;close all;clear all;

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

type = 3;
n = 300;

%% input test
data = dataCase(type, n);
figure
plot(data(1,1:n),data(2,1:n),'ro');hold on;
plot(data(1,n+1:2*n),data(2,n+1:2*n),'bo');

labels{1} = [1.*ones(1,n),2.*ones(1,n)];
labels{2} = [1.*ones(1,n),2.*ones(1,n)];
labels{3} = [1.*ones(1,2*n),2.*ones(1,2*n)];
label = labels{type};

% data = dataCase(2, n);
% figure
% plot(data(1,1:n),data(2,1:n),'r.');hold on;
% plot(data(1,n+1:2*n),data(2,n+1:2*n),'b.');
% 
% data = dataCase(3, n);
% figure
% plot(data(1,1:2*n),data(2,1:2*n),'r.');hold on;
% plot(data(1,2*n+1:4*n),data(2,2*n+1:4*n),'b.');


%% Q1 test with random weights
num_of_hidden_units = [10 20 10];% mxn: m is the hidden units per layer, n - layer
num_inputs = size(data,1);
num_outputs = 2;

htb = zeros(1,size(num_of_hidden_units,1)+1);
ztb = zeros(1,size(num_of_hidden_units,1)+1);

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
    
    htb(i) = column-1;
    ztb(i) = row;
end

x = data;
x = normalizeData(x);
[y,h,d] = forwardPropagation(nn,x,@ReLU);
[val,id]=max(y);
figure
subplot(1,2,1)
plot(data(1,label == 1),data(2,label == 1),'r.');hold on;
plot(data(1,label == 2),data(2,label == 2),'b.');
title('Truth');
set(gca,'FontName','Arial','FontSize',15);
subplot(1,2,2)
plot(data(1,id == 1),data(2,id == 1),'r.');hold on;
plot(data(1,id == 2),data(2,id == 2),'b.');
title('Result: Random Wight')
set(gca,'FontName','Arial','FontSize',15);

%% training using SGD
oldloss = -1e-6;
epsilon1 = 1e-6;

cyclic_id = 1;

datalength = size(x,2);
iter = 0;
losses = [];

t = [ones(1,size(x,2));2.*ones(1,size(x,2))];
id = (t(1,:)==label);t(1,id) = 1; t(1,~id) = 0; 
t(2,:) = 1 - t(1,:);

while iter < 1e4
    % forward
    [y,h,z] = forwardPropagation(nn,x,@ReLU);
    % compute loss
    newloss = loss(t,y);
    if abs(newloss - oldloss) < 1e-10
        break;
    end
    disp(newloss);
    oldloss = newloss;
    % SGD: backpropagation
    hc = cell2mat(h);hc = hc(:,cyclic_id);hc = mat2cell(hc,htb);
    zc = cell2mat(z);zc = zc(:,cyclic_id);zc = mat2cell(zc,ztb);
    grad = backpropagation(nn,t(:,cyclic_id),y(:,cyclic_id),hc,zc,@ReLUDer);
    % SDG: gradient descent
    nn = cellfun(@updateW,nn,grad,'UniformOutput', false);
    
%     cyclic_id = cyclic_id + 1;if cyclic_id > datalength; cyclic_id = 1;end
    cyclic_id = randperm(datalength,1);
    
    iter = iter + 1;
    losses(iter) = newloss;
end

figure
plot(losses);
grid on;
title('Loss');

[y,h,z] = forwardPropagation(nn,x,@ReLU);
[val,id]=max(y);
figure
subplot(1,2,1)
plot(data(1,label == 1),data(2,label == 1),'r.');hold on;
plot(data(1,label == 2),data(2,label == 2),'b.');
title('Truth');
set(gca,'FontName','Arial','FontSize',15);
subplot(1,2,2)
plot(data(1,id == 1),data(2,id == 1),'r.');hold on;
plot(data(1,id == 2),data(2,id == 2),'b.');
title('Result: Random Wight')
set(gca,'FontName','Arial','FontSize',15);

function W = updateW(W,grad)
    learning_rate = 0.01;
    W = W-grad.*learning_rate;
end

function data = dataCase(type, n)
% setting up data, 3 types. n is the number of points for each set.
    if type == 1
        r = rand([1,n]).*5;
        theta = linspace(0,2*pi,n);
        pt = [r.*cos(theta);r.*sin(theta)];
        d = 6;
        sx = d*cos(pi/4);sy = d*sin(pi/4);
        data = [pt+repmat([sx;sy],1,n) pt-repmat([sx;sy],1,n)];
    elseif type == 2
        r1 = rand([1,n]).*5;
        theta = linspace(0,2*pi,n);
        pt1 = [r1.*cos(theta);r1.*sin(theta)];
        r2 = r1 + 6;
        pt2 = [r2.*cos(theta);r2.*sin(theta)];
        data = [pt2 pt1];
    else
        r = rand([1,n]).*5;
        theta = linspace(0,2*pi,n);
        pt = [r.*cos(theta);r.*sin(theta)];
        d = 7;
        sx1 = d*cos(pi/4);sy1 = d*sin(pi/4);
        sx2 = d*cos(3*pi/4);sy2 = d*sin(3*pi/4);
        data = [pt+repmat([sx1;sy1],1,n) pt-repmat([sx1;sy1],1,n) pt+repmat([sx2;sy2],1,n) pt-repmat([sx2;sy2],1,n)];
    end
end


function dn = normalizeData(d)
% d is assume to be mxn, m is the dimention, n is the number
    dm = mean(d,2);
    dc = d - repmat(dm,1,size(d,2));
    scale = sqrt(1 ./ mean(diag(dc'*dc)));
    dn = dc .* scale;
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
    A1 = [h{n}' ones(size(y,2),1)];
    A2 = delta{n};
    grad{n} = A1.*A2;
    for i = length(W)-1:-1:1
        d1 = W{i+1}'*delta{i+1};
        delta{i} = fgrad(z{i}).*d1(1:end-1);
        grad{i} = [h{i}' 1].*delta{i};
    end
end

%% todo minibatch. Plan: 1. SGD; 2. Batch; 3. Minibatch