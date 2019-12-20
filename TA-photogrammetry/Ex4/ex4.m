clc;close all;clear all;

%% Harris corner detector
imColor=imread('./data/library2.jpg');
figure
imshow(imColor)
if size(imColor,3) == 3
    imGray=rgb2gray(imColor);
else
    imGray = imColor;
end
im = im2double(imGray);

% gradient
Ix = edge(im,'Sobel',[],'horizontal');
Iy = edge(im,'Sobel',[],'vertical');
figure;
subplot(1,2,1);imshow(Ix);
subplot(1,2,2);imshow(Iy);

Ixx = Ix.*Ix;
Iyy = Iy.*Iy;
Ixy = Ix.*Iy;
g = fspecial('gaussian', 3, 1);% try play with this two values to see what will hapen
Ixx = imfilter(Ixx,g,'replicate');
Iyy = imfilter(Iyy,g,'replicate');
Ixy = imfilter(Ixy,g,'replicate');
figure;
subplot(1,3,1);imshow(Ixx);
subplot(1,3,2);imshow(Iyy);
subplot(1,3,3);imshow(Ixy);

C = Ixx.*Iyy - Ixy.^2 - 0.04.*(Ixx+Iyy).^2;

threshold = 0.3*max(C(:));
C(C(:)<threshold) = 0;
figure;
imshow(C);
 
[row,col] = nonmaxsuppts(C,'radius', 2, 'N', 1000);

figure
img=imshow(imColor),title('my-Harris'),
hold on
plot(col,row, 'ro','MarkerSize',10),
hold off
%% this ends the Harris part

%% RANSAC part
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






