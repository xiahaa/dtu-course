clc;close all;clear all;
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end
%% Q1 
qhomo = [[3 4 1]; ...
     [4 6.5 0.5]; ...
     [5 2 0.1]; ...
     [10 100 20]; ...
     [0.1 0.2 100]]';
qinhomo = qhomo ./ qhomo(end,:);

%% Q2
qhomo = [[3 4 1 0.1]; ...
     [4 6.5 0.5 2]; ...
     [5 2 0.1 3]; ...
     [10 100 20 1]; ...
     [0.1 0.2 100 10]]';
qinhomo = qhomo ./ qhomo(end,:);


%% Q3
p1 = [1 2 1]';p2 = [5 3.5 1]';
l1 = cross(p1,p2);scale = 1./norm(l1(1:2));
l1 = l1.*scale;

%% Q4
p3 = [7.5 3.6 1]';
dist = abs(p3'*l1);


%% Q5
l1 = [1 1 -1]';
l2 = [-1 3 -4]';
pinter = cross(l1,l2);
pinter = pinter ./ pinter(3);

%% Q6: Euclidean transformation 2D, length, angle, orthogonal, parallel are preserved
x = 0:0.1:5;
y = 0:0.1:5;
pt = [x repmat(x(end), 1, numel(y)) fliplr(x) repmat(x(1), 1, numel(y)); ...
      repmat(y(1), 1, numel(y)) y repmat(y(end), 1, numel(y)) fliplr(y)];
pt = [pt;ones(1,size(pt,2))];

figure
plot(pt(1,:), pt(2,:), 'LineWidth', 2);hold on;
a = 30*pi/180.0;
t = [2;-2];
A = [cos(a)    sin(a)         t(1); ...
    -sin(a)    cos(a)          t(2); ...
      0 0      1.0000];
pt1 = A*pt;
pt1 = pt1./pt1(3,:);
plot(pt1(1,:), pt1(2,:),'r--');hold on;axis equal;grid on;

%% Q7: similarity transformation, angle, orthogonal, parallel, ratio of length are preserved
close all;
s = 0.2;
A = [cos(a)*s    s*sin(a)         t(1); ...
    -sin(a)*s    s*cos(a)         t(2); ...
      0 0      1.0000];
pt1 = A*pt;
pt1 = pt1./pt1(3,:);
plot(pt(1,:), pt(2,:), 'LineWidth', 2);hold on;
plot(pt1(1,:), pt1(2,:),'r--');hold on;axis equal;grid on;
    
  
%% Q8: affine 2D, parallel, ratio of area, see shearing effect
close all;
A = [0.5 0.3 2; ...
     0.1 1.2 -2; ...
     0 0 1];  
pt1 = A*pt;
pt1 = pt1./pt1(3,:);
plot(pt(1,:), pt(2,:), 'LineWidth', 2);hold on;
plot(pt1(1,:), pt1(2,:),'r--');hold on;

%% Q9:
clc;clear all;close all;
Q=Box3D;
plot3(Q(1,:),Q(2,:),Q(3,:),'.');
axis equal;
axis([-1 1 -1 1 -1 5]);
xlabel('x');
ylabel('y');
zlabel('z');
K = [1 0 0;0 1 0;0 0 1];
Qh = Q ./ Q(3,:);
q = K*Qh;
plot(q(1,:),q(2,:),'.')
axis equal
% axis([0 640 0 480])
a11 = 0.2;a2 = 0.1;a3 = 0.0;a4 = 0;
r = (q(1,:).^2+q(2,:).^2);
fr = 1+r.*a11+r.^2.*a2+r.^3.*a3+r.^4.*a4;
q(1:2,:) = q(1:2,:).*fr;
% q = K*q;
% q = q./q(3,:);
hold on;
plot(q(1,:),q(2,:),'.')
axis equal

%% Q10: 
I = imread('Tiles_perspective_distort.png');
I = im2double(I);
% I = rgb2gray(I);
h = figure;
imshow(I);
[x,y] = ginput(4);
q = round([x y]');
q = [q;ones(1,4)];

l1 = cross(q(:,1),q(:,2));
l2 = cross(q(:,3),q(:,4));
l3 = cross(q(:,1),q(:,4));
l4 = cross(q(:,2),q(:,3));

I = drawlines(I,q,[[1 2];[3 4];[1 4];[2 3]]);
imshow(I);

%
p1 = cross(l1,l2);p1 = p1./p1(3);
p2 = cross(l3,l4);p2 = p2./p2(3);

%
linf = cross(p1,p2);
linf = linf ./ linf(3);

H = [1 0 0;0 1 0;linf(1) linf(2) linf(3)];

%% try play with other transformation
% H = [1.1 -1.1 0;1.1 1.1 0;0 0 1];
%% e
% H = [cos(-5/57.3) sin(-5/57.3) 40;-sin(-5/57.3) cos(-5/57.3) 100;0 0 1];
%% s
% s = 1
% H = [s*cos(-5/57.3) s*sin(-5/57.3) 0;-s*sin(-5/57.3) s*cos(-5/57.3) 0;0 0 1];
%% affine
% H = [s*cos(-5/57.3) s*sin(-5/57.3) 2;-s*sin(-5/57.3) s*cos(-5/57.3) 2;0 0 1];
Iw = warpping(I,H);

if size(Iw,3) == 3
    Iwr = Iw(:,:,1);Iwr = imresize(Iwr,[size(I,1),size(I,2)]);
    Iwg = Iw(:,:,2);Iwg = imresize(Iwg,[size(I,1),size(I,2)]);
    Iwb = Iw(:,:,3);Iwb = imresize(Iwb,[size(I,1),size(I,2)]);
    Iw = cat(3,Iwr,Iwg,Iwb);
else
    Iw = imresize(Iw,size(I));
end

Ig = cat(2,I,Iw);
imshow(Ig)