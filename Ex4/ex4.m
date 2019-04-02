clc;close all;clear all;

[X, lineTrue] = gen_line_data(500);

ph = tohomo(X);

ransac.pinlier = 0.99;
ransac.estt_fun = @line_estimation;%plane_estimation
ransac.eval_fun = @dist2line;%dist2plane
ransac.maxiter = 1e6;
ransac.threshold = 0.1;
ransac.inliers = [];
ransac.minimumset = 2;
result = ransac_routine(ph, ransac);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Results Visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on
% ind = results.CS;
plot(X(1, result.inliers), X(2, result.inliers), 'sg','MarkerSize', 5);
plot(X(1, ~result.inliers), X(2, ~result.inliers), 'sr','MarkerSize', 5);
xmin = min(X(1,:));
xmax = max(X(1,:));
xx = linspace(xmin,xmax,100);
yytrue = -(lineTrue(1).*xx+lineTrue(3))./(lineTrue(2)+1e-6);
yy = -(result.params(1).*xx+result.params(3))./(result.params(2)+1e-6);
plot(xx, yytrue, 'k-.','LineWidth',1);
plot(xx, yy, 'm--','LineWidth',2);
legend('Estimated Iniliers', 'Estimated Outliers','True Line','Estimated Line');
xlabel('x')
ylabel('y')
title('RANSAC results for 2D line estimation')
axis equal tight
set(gca,'FontName','Arial','FontSize',20);






