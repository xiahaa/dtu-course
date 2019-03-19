clc;close all;clear all;
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

%% Q1-2:
clc;clear all;close all;
Q=Box3D;
plot3(Q(1,:),Q(2,:),Q(3,:),'.');
axis equal;
axis([-1 1 -1 1 -1 5]);
xlabel('x');
ylabel('y');
zlabel('z');

Qh = tohomogeneous(Q);

R = Rxyz(0.2, -0.3, 0.1);
t = [0.88;0.57;0.19];
f = 1000; cu = 300; cv = 200;
K = [f 0 cu;0 f cv;0 0 1];
P = K*[R t];
qc = P*Qh;
q = qc ./ qc(3,:);
plot(q(1,:),q(2,:),'.')
axis equal
axis([0 640 0 480])

k1 = -5e-1;k2 = -3e-1;k3 = -5e-1;
p1 = 0;%3e-2;
p2 = 0;%-2e-2; 

qn = K\q;

x = qn(1,:); y = qn(2,:);
r = (x.^2+y.^2);
dradial = 1+r.*k1+r.^2.*k2+r.^3.*k3;
dtangentx = 2*p1.*x.*y + p2.*(r + 2.*x.^2);
dtangenty = p1.*(r + 2.*y.^2) + 2*p2.*x.*y;

qd(1,:) = x.*dradial + dtangentx;
qd(2,:) = y.*dradial + dtangenty;
qd(3,:) = ones(1,size(qn,2));
q1 = K*qd;
q1 = q1./q1(3,:);
hold on;
plot(q1(1,:),q1(2,:),'.')
axis equal

% Q3
% Create a set of calibration images.
% images = imageDatastore('./calibration');
% imageFileNames = images.Files;
% 
% % Detect calibration pattern.
% [imagePoints, boardSize] = detectCheckerboardPoints(imageFileNames);
% 
% % Generate world coordinates of the corners of the squares.
% squareSize = 112; % millimeters
% worldPoints = generateCheckerboardPoints(boardSize, squareSize);
% 
% % Calibrate the camera.
% I = readimage(images, 1); 
% imageSize = [size(I, 1), size(I, 2)];
% [params, ~, estimationErrors] = estimateCameraParameters(imagePoints, worldPoints, ...
%                                      'ImageSize', imageSize);

load('gopro_27.mat');
I = imread('distort.bmp');
K = cameraParams.IntrinsicMatrix';
k1 = cameraParams.RadialDistortion(1);
k2 = cameraParams.RadialDistortion(2);
k3 = cameraParams.RadialDistortion(3);
p1 = cameraParams.TangentialDistortion(1);
p2 = cameraParams.TangentialDistortion(2);
Irec = undistortImg(I, K, k1, k2, k3, p1, p2);

% resize and show
im1 = imresize(I,0.4);
im2 = imresize(Irec,0.4);
imc = cat(2,im1,im2);
figure
imshow(imc);

% test undistort point
figure
imshow(I);
[x,y] = ginput(10);
x1 = round([x y]');
% plot(x1(1,:),x1(2,:),'r-','MarkerSize',10);
Ic = cat(3,I,I,I);
indices = zeros(size(x1,2)-1,2);
for i = 1:size(x1,2)-1
    indices(i,:) = [i,i+1];
end
for i = 1:size(indices,1)
    Ic = drawline(Ic,[x1(1,indices(i,1)) x1(1,indices(i,2))], [x1(2,indices(i,1)) x1(2,indices(i,2))], [255,0,0]);
end

x2 = undistortPoint(x1, K, k1, k2, k3, p1, p2);
Irecc = cat(3,Irec,Irec,Irec);
for i = 1:size(indices,1)
    Irecc = drawline(Irecc,[x2(1,indices(i,1)) x2(1,indices(i,2))], [x2(2,indices(i,1)) x2(2,indices(i,2))], [0,255,0]);
end

Ishow = cat(2,Ic,Irecc);
imshow(Ishow);


% I = drawlines(I,q,[[1 2];[3 4];[1 4];[2 3]]);







