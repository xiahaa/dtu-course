clc;clear;close all;
%% exercise 5 of bayesian scientifical computing
xreal0 = [61000, 3048, 19161]';
j = 0:1:300;
deltat = 0.1;
t = deltat.*j;

[t,x] = ode45(@func_ode, t, xreal0);
x = x';
q = [1 0 0]';

sigma = 100;
colormap = 'jet';

b = q'*x;
bn = b + sigma.*randn(size(b));
figure(1);
plot(t, b, 'LineWidth', 2);
hold on; grid on;
plot(t, bn, 'LineWidth', 2);
legend({'true','noise'},'FontSize',15);

hinit = b(1);
vinit = 3000;
betainit = 20000;

%% initialization
x0 = [hinit, vinit, betainit]';
N = 5000;
C = diag([500,200,1500]);
C1 = diag([100,100,5]);
D = diag([100,100,500]);
xp = x0 + C * randn(3,N);
wp = ones(1,N)./N;

xpost = zeros(length(x0),length(t));
xmedian = zeros(length(x0),length(t));
x_25 = zeros(length(x0),length(t));
x_75 = zeros(length(x0),length(t));
xpost(:,1) =  sum(xp.*wp,2);
xmedian(:,1) = median(xp,2);
[x_1,x_2] = func_p_rule(xp, 25);
x_25(:,1) = x_1;
x_75(:,1) = x_2;

stdpost = zeros(length(x0), length(t));
stdpost(:,1) = diag(C);

dt = deltat;
figure
for i = 1:1:length(t)
    %% likelihood
    lh = exp(-0.5/sigma^2.*(xp(1,:)-bn(i)).^2);
    wpn = wp.*lh;
    wp = wpn / sum(wpn);
    cumwp = cumsum(wp);
    %% random draw with replacement, one way
    pdraw = rand(1,N);
    indices = arrayfun(@(pd) (func_min(pd,cumwp)), pdraw);
    plot(sort(indices));pause(0.1);
    xp = xp(:,indices);
    xpnew = xp + D*randn(3,N);
    %% update weights
    wp = exp(-0.5/sigma^2.*((xpnew(1,:)-bn(i)).^2 - (xp(1,:)-bn(i)).^2));
    wp = wp ./ sum(wp);
    xp = xpnew;
    
    %% posterior
    xpost(:,i) = sum(xp.*wp,2);
    err = xp - xpost(:,i);%3xN
    Cov = (err.*wp)*err';
    stdpost(:,i) = sqrt(diag(Cov));
    xmedian(:,i) = median(xp,2);
    [x_1,x_2] = func_p_rule(xp, 25);
    x_25(:,i) = x_1;
    x_75(:,i) = x_2;
    
    dx = func_ode(t, xp);
    %% prior
    xp = xp + dt.*dx;% + C1 * randn(3,N);
end

figure(2);
plot(t,x(1,:),'b-','LineWidth',2);hold on;grid on;
plot(t,xpost(1,:),'r-','LineWidth',2);
plot(t,xmedian(1,:),'g-','LineWidth',2);
plot(t,x_25(1,:),'m--','LineWidth',2);
plot(t,x_75(1,:),'m--','LineWidth',2);
legend({'ode','mean','median','25%','75%'},'FontSize',15);

figure(3);
plot(t,x(2,:),'b-','LineWidth',2);hold on;grid on;
plot(t,xpost(2,:),'r-','LineWidth',2);
plot(t,xmedian(2,:),'g-','LineWidth',2);
plot(t,x_25(2,:),'m--','LineWidth',2);
plot(t,x_75(2,:),'m--','LineWidth',2);
legend({'ode','mean','median','25%','75%'},'FontSize',15);

figure(4);
plot(t,x(3,:),'b-','LineWidth',2);hold on;grid on;
plot(t,xpost(3,:),'r-','LineWidth',2);
plot(t,xmedian(3,:),'g-','LineWidth',2);
plot(t,x_25(3,:),'m--','LineWidth',2);
plot(t,x_75(3,:),'m--','LineWidth',2);
plot(t,xreal0(3)*ones(1,length(t)),'k-','LineWidth',2);
legend({'ode','mean','median','25%','75%','real'},'FontSize',15);

function [x_1,x_2] = func_p_rule(x, p)
    for i = 1:size(x,1)
        [val,id] = sort(x(i,:),'ascend');
        x_1 = x(:, round((1-p/100)*length(x)*0.5));
        x_2 = x(:, round((1+p/100)*length(x)*0.5));
    end
end

function id = func_min(x, wp)
    id = find((wp-x)>=0,1);
end

function dx = func_ode(t, x)
    g = 9.81;%m/s^2
    gamma = 1.754;%kg/m^3
    eta = 1.39*1e-4;%m^-1
    rho = @(h) (gamma.*exp(-eta.*h));
    dx = x;
    dx(1,:) = -x(2,:);
    dx(2,:) = g - g.*rho(x(1,:)).*(x(2,:).^2)./(2*x(3,:));
    dx(3,:) = 0;
end