clc;close all;clear all;

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

type = 1;
n = 1000;

%% input test
data = dataCase(type, n);
labels{1} = [1.*ones(1,n),2.*ones(1,n)];
labels{2} = [1.*ones(1,n),2.*ones(1,n)];
labels{3} = [1.*ones(1,2*n),2.*ones(1,2*n)];
label = labels{type};

figure
id = label == 1;
plot(data(1,id),data(2,id),'ro');hold on;
plot(data(1,~id),data(2,~id),'bo');

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
num_of_hidden_units = [10 10];% mxn: m is the hidden units per layer, n - layer
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
[x,xm,xscale] = normalizeData(x);
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

datalength = size(x,2);
iter = 0;
losses = [];

t = [ones(1,size(x,2));2.*ones(1,size(x,2))];
id = (t(1,:)==label);t(1,id) = 1; t(1,~id) = 0; 
t(2,:) = 1 - t(1,:);

epoches = 200;
batchSize = 100;
learning_rate = 0.01;
old_grad = {};
momentum = 0.9;

dummy1 = 1:(batchSize):size(x,2);
if dummy1(end) ~= size(x,2)
    dummy1(end+1) = size(x,2);
end
batcheIndex = [[dummy1(1:end-1)]' [dummy1(2:end-1)'-1;dummy1(end)]];
figure;grid on;hold on;
tic
for i = 1:epoches
    rd = randperm(size(x,2),size(x,2));
    for j = 1:size(batcheIndex,1)
        id = rd(batcheIndex(j,1):batcheIndex(j,2));
        xb = x(:,id);
        tb = t(:,id);
        % forward
        [y,h,z] = forwardPropagation(nn,xb,@ReLU);
        % SGD: backpropagation
        grad = backpropagation(nn,tb,y,h,z,@ReLUDer);
        % SDG: gradient descent
%       nn = cellfun(@updateW,nn,grad,'UniformOutput', false);
        [nn, old_grad] = updateSGD(nn, grad, learning_rate, momentum, old_grad);
    end
    % compute loss
    [y,~,~] = forwardPropagation(nn,x,@ReLU);
    newloss = loss(t,y);
    if abs(newloss - oldloss) < 1e-6
        break;
    end
    %disp(newloss);
    oldloss = newloss;
    losses(i) = newloss;
end
toc
plot(losses);
title('Loss');
xlabel('Epoches');
ylabel('Loss');

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

%% for all region
xmin = floor(min(data(1,:)));
xmax = round(max(data(1,:)));
ymin = floor(min(data(2,:)));
ymax = round(max(data(2,:)));

xx = xmin:0.1:xmax;
yy = ymin:0.1:ymax;

[xg,yg] = meshgrid(xx,yy);
xg = xg(:);
yg = yg(:);
xgg = cat(1,xg',yg');
xgg = xgg - repmat(xm,1,size(xgg,2));
xgg = xgg .* xscale;

[ygg,~,~] = forwardPropagation(nn,xgg,@ReLU);
[val,id]=max(ygg);
im = zeros(length(xx),length(yy));
im(id == 1) = 0.4;
im(id == 2) = 0.8;
imc = cat(3,im,im,im);

data1 = data(:,label==1);
xdiff = abs(repmat(data1(1,:),length(xx),1)-repmat(xx',1,size(data1,2)));
[~,id1] = min(xdiff);
ydiff = abs(repmat(data1(2,:),length(yy),1)-repmat(yy',1,size(data1,2)));
[~,id2] = min(ydiff);
id = (id2-1).*length(xx)+id1;
imr = imc(:,:,1);imr(id) = 1;
img = imc(:,:,2);img(id) = 0;
imb = imc(:,:,3);imb(id) = 0;
imc = cat(3,imr,img,imb);

data2 = data(:,label==2);
xdiff = abs(repmat(data2(1,:),length(xx),1)-repmat(xx',1,size(data2,2)));
[~,id1] = min(xdiff);
ydiff = abs(repmat(data2(2,:),length(yy),1)-repmat(yy',1,size(data2,2)));
[~,id2] = min(ydiff);
id = (id2-1).*length(xx)+id1;
imr = imc(:,:,1);imr(id) = 0;
img = imc(:,:,2);img(id) = 0;
imb = imc(:,:,3);imb(id) = 1;

imc = cat(3,imr,img,imb);
imc = flipud(imc);
figure;imshow(imc);


function [W,old_grad] = updateSGD(W, grad, learning_rate, momentum, old_grad)
    if isempty(old_grad)
        for i = 1:length(grad)
            old_grad{i} = -grad{i}.*learning_rate;
        end
    else
        for i = 1:length(grad)
            old_grad{i} = old_grad{i}.*momentum - grad{i}.*learning_rate;
        end
    end
    for i = 1:length(W)
        W{i} = W{i}+old_grad{i};
    end
end

function W = updateW(W,grad)
    learning_rate = 0.1;
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


function [dn,dm,scale] = normalizeData(d)
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
    w = randn([m,1]).*sqrt(2/n);
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