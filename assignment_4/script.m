clc;close all;clear all;

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

type = 1;
n = 500;

%% input test
data = dataCase(1, n);
figure
plot(data(1,1:n),data(2,1:n),'r.');hold on;
plot(data(1,n+1:2*n),data(2,n+1:2*n),'b.');

t = [1.*ones(1,n),2.*ones(1,n)];

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
num_of_hidden_units = [4];% mxn: m is the hidden units per layer, n - layer
num_inputs = size(data,1);
num_outputs = 2;

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
end

x = data;
x = normalizeData(x);
y = forwardPropagation(nn,x,@ReLU);
[val,id]=max(y);
figure
subplot(1,2,1)
plot(data(1,t == 1),data(2,t == 1),'r.');hold on;
plot(data(1,t == 2),data(2,t == 2),'b.');
title('Truth');
set(gca,'FontName','Arial','FontSize',15);
subplot(1,2,2)
plot(data(1,id == 1),data(2,id == 1),'r.');hold on;
plot(data(1,id == 2),data(2,id == 2),'b.');
title('Result: Random Wight')
set(gca,'FontName','Arial','FontSize',15);


function data = dataCase(type, n)
% setting up data, 3 types. n is the number of points for each set.
    if type == 1
        r = rand([1,n]).*5;
        theta = linspace(0,2*pi,n);
        pt = [r.*cos(theta);r.*sin(theta)];
        d = 5;
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
        d = 6.5;
        sx1 = d*cos(pi/4);sy1 = d*sin(pi/4);
        sx2 = d*cos(3*pi/4);sy2 = d*sin(3*pi/4);
        data = [pt+repmat([sx1;sy1],1,n) pt-repmat([sx1;sy1],1,n) pt+repmat([sx2;sy2],1,n) pt-repmat([sx2;sy2],1,n)];
    end
end


function dn = normalizeData(d)
% d is assume to be mxn, m is the dimention, n is the number
    dm = mean(d,2);
    dc = d - repmat(dm,1,size(d,2));
    scale = sqrt(2 ./ mean(diag(dc'*dc)));
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
    hdz(z<=0) = 0;
end

function L = loss(t,y)
% summation of summation along each class. 
    L = -sum(t.*log(y));
end

function [y,h,z] = forwardPropagation(W,x,f)
% run forward propagation for the forward nn. W is assumed to be a cell
% with each element being the weight matrix of layer i. f is assumed to be
% the activation function handler.
    z = cell(length(W),1);
    h = cell(length(W),1);
    
    h{1} = [x;ones(1,size(x,2))];
    
    z = W{1}*;
    z = f(z);
    for i = 1:length(W)-1
        % layer to layer
        z{i} = W{i}*[h{i};1];
        h{i+1} = f(z{i});
        
        znew = W{i}*[z;ones(1,size(z,2))];
        z = f(znew);
    end
    y = softMax(z);
end

function grad = backpropagation(W,t,y,fgrad)
% run backpropagation for gradient computation. W is assumed to be a cell.
    % first, take care of output
    grad = cell(length(W),1);
    delta =  y - t;
    grad
    
end

%% todo minibatch. Plan: 1. SGD; 2. Batch; 3. Minibatch