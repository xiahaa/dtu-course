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

%% Q3-Q4
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
% figure
% imshow(I);
% [x,y] = ginput(10);
% x1 = round([x y]');
% % plot(x1(1,:),x1(2,:),'r-','MarkerSize',10);
% Ic = cat(3,I,I,I);
% indices = zeros(size(x1,2)-1,2);
% for i = 1:size(x1,2)-1
%     indices(i,:) = [i,i+1];
% end
% for i = 1:size(indices,1)
%     Ic = drawline(Ic,[x1(1,indices(i,1)) x1(1,indices(i,2))], [x1(2,indices(i,1)) x1(2,indices(i,2))], [255,0,0]);
% end
% 
% x2 = undistortPoint(x1, K, k1, k2, k3, p1, p2);
% Irecc = cat(3,Irec,Irec,Irec);
% for i = 1:size(indices,1)
%     Irecc = drawline(Irecc,[x2(1,indices(i,1)) x2(1,indices(i,2))], [x2(2,indices(i,1)) x2(2,indices(i,2))], [0,255,0]);
% end
% 
% Ishow = cat(2,Ic,Irecc);
% imshow(Ishow);

% % I = drawlines(I,q,[[1 2];[3 4];[1 4];[2 3]]);


%% Q5-Q6
% N = 100;
% p = rand([3,N]) * 5 - 2.5;
% 
% C1 = [10;0;0];
% C2 = [0;10;0];
% C3 = [0;0;10];
% 
% R1 = [0 1 0;0 0 -1;-1 0 0];
% R2 = [-1 0 0;0 0 -1;0 -1 0];
% R3 = [0 1 0;1 0 0;0 0 -1];
% t1 = -R1*C1;
% t2 = -R2*C2;
% t3 = -R3*C3;
% 
% figure
% plot3(p(1,:),p(2,:),p(3,:),'g*');hold on;
% cam1 = plotCamera('Location',C1,'Orientation',R1,'Opacity',0,'Color',[1 0 0]);
% cam2 = plotCamera('Location',C2,'Orientation',R2,'Opacity',0,'Color',[0 1 0]);
% cam3 = plotCamera('Location',C3,'Orientation',R3,'Opacity',0,'Color',[0 0 1]);
% xlabel('x');
% ylabel('y');
% zlabel('z');
% 
% K = [250 0 160;0 250 120;0 0 1];
%     
% [uv1, in1] = proj(R1,t1, p, K);
% [uv2, in2] = proj(R2,t2, p, K);
% [uv3, in3] = proj(R3,t3, p, K);
% 
% im = zeros(240,320);
% figure;imshow(im);hold on;
% plot(uv1(2,in1),uv1(1,in1),'ro');
% figure;imshow(im);hold on;
% plot(uv2(2,in2),uv2(1,in2),'go');
% figure;imshow(im);hold on;
% plot(uv3(2,in3),uv3(1,in3),'bo');
% 
% % for i = 
% q1 = uv1(:,in1);
% p1 = p(:,in1);
% 
% q2 = uv2(:,in2);
% p2 = p(:,in2);
% 
% q3 = uv3(:,in3);
% p3 = p(:,in3);
clear all;
load('p3p.mat');

[Rc1,tc1] = p3p_Grunert(p1(:,1:3), q1(:,1:3), K);
[Rc2,tc2] = p3p_Grunert(p2(:,1:3), q2(:,1:3), K);
[Rc3,tc3] = p3p_Grunert(p3(:,1:3), q3(:,1:3), K);
    
[Rc1,tc1] = selectBestPose(Rc1,tc1,p1,q1,K);
[Rc2,tc2] = selectBestPose(Rc2,tc2,p2,q2,K);
[Rc3,tc3] = selectBestPose(Rc3,tc3,p3,q3,K);

Cc1 = -Rc1'*tc1;
Cc2 = -Rc2'*tc2;
Cc3 = -Rc3'*tc3;

C1 = -R1'*t1;
C2 = -R2'*t2;
C3 = -R3'*t3;

figure
% plot3(p1(1,:),p(2,:),p(3,:),'g*');
cam1 = plotCamera('Location',C1,'Orientation',R1,'Opacity',0.0,'Color',[1 0 0],'Label','Camera1');hold on;
cam2 = plotCamera('Location',C2,'Orientation',R2,'Opacity',0.0,'Color',[1 0 0],'Label','Camera2');
cam3 = plotCamera('Location',C3,'Orientation',R3,'Opacity',0.0,'Color',[1 0 0],'Label','Camera3');

camc1 = plotCamera('Location',Cc1,'Orientation',Rc1,'Opacity',0.2,'Color',[0 1 0],'Size',0.5,'Label','Est1');
camc2 = plotCamera('Location',Cc2,'Orientation',Rc2,'Opacity',0.2,'Color',[0 1 0],'Size',0.5,'Label','Est2');
camc3 = plotCamera('Location',Cc3,'Orientation',Rc3,'Opacity',0.2,'Color',[0 1 0],'Size',0.5,'Label','Est3');
xlabel('x:(m)','FontName','Aerial','FontSize',15);
ylabel('y:(m)','FontName','Aerial','FontSize',15);
zlabel('z:(m)','FontName','Aerial','FontSize',15);

function [Ropt,topt] = selectBestPose(R,t,p,q,K)
    minerr = 1e6;
    minid = 0;
    ph = [p;ones(1,size(p,2))];
    for i = 1:size(R,3)
        P1 = K*([R(1:3,1:3,i) t(1:3,1,i)]);
        uv1rep = P1*ph;
        uv1rep = uv1rep./uv1rep(3,:);
        err = uv1rep(1:2,:) - q(1:2,:);
        avgerr = sum(diag(err'*err)) / size(q,2);
        if avgerr < minerr
            minerr = avgerr;
            minid = i;
        end
    end
    Ropt = R(:,:,minid);
    topt = t(:,:,minid);
end





