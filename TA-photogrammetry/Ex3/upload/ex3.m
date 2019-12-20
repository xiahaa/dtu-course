%% loading-8pointalgorithm Q1-Q2
im1 = rgb2gray(imread('000002.bmp'));
im2 = rgb2gray(imread('000003.bmp'));
load('calib.mat');
[im1,~] = undistortImage(im1,cameraParams);
[im2,~] = undistortImage(im2,cameraParams);

load('Fdata.mat');

im3 = cat(2,im1,im2);
figure;imshow(im3);hold on;

plot(x1(:,1),x1(:,2),'ro');
plot(x2(:,1)+size(im1,2),x2(:,2),'go');

shift = size(im1,2);
cmap = lines(5);
k = 1;
for i = 1:size(x1,1)
    ptdraw = [x1(i,1), x1(i,2);
              x2(i,1)+shift, x2(i,2)];
    plot(ptdraw(:,1),ptdraw(:,2),'LineStyle','-','LineWidth',1,'Color',cmap(k,:));
    k = mod(k+1,5);if k == 0 k = 1;end
end
A=[ptdraw(1,1)*ptdraw(1,2) ptdraw(1,1)*ptdraw(2,2) ptdraw(1,1) ptdraw(2,1)*ptdraw(1,2) ptdraw(2,1)*ptdraw(2,2) ptdraw(2,1) ptdraw(1,2) ptdraw(2,2) 1];
[U S V]=svd(A);
F = V(:,end);
vgg_gui_F(im1,im2,F');
%triangulation Q3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Q3   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N = 100;
% p = rand([3,N]) * 10 - 5;% from near to far
% p(1,:) = p(1,:) + 15;
% C1 = [0;0;0];
% C2 = [0;0.2;0];
% 
% rad1 = -10/57.3;
% R1 = [cos(rad1) 0 -sin(rad1);0 1 0;sin(rad1) 0 cos(rad1)]*[0 -1 0;0 0 -1;1 0 0];
% R2 = [cos(-rad1) 0 -sin(-rad1);0 1 0;sin(-rad1) 0 cos(-rad1)]*[0 -1 0;0 0 -1;1 0 0];
% 
% t1 = -R1*C1;
% t2 = -R2*C2;
% 
% figure
% h1 = plot3(p(1,:),p(2,:),p(3,:),'g*');hold on;
% cam1 = plotCamera('Location',C1,'Orientation',R1,'Opacity',0,'Color',[1 0 0],'Size',0.4,'Label','Camera1');
% cam2 = plotCamera('Location',C2,'Orientation',R2,'Opacity',0,'Color',[0 1 0],'Size',0.4,'Label','Camera2');
% axis equal
% xlabel('x:(m)');
% ylabel('y:(m)');
% zlabel('z:(m)');
% title('Triangulation Simulation');
% set(gca,'FontName','Arial','FontSize',20);
% 
% K = [500 0 320;0 500 240;0 0 1];
% 
% [uv1, in1] = proj(R1,t1, p, K);
% [uv2, in2] = proj(R2,t2, p, K);
% 
% in = in1 & in2;
% q1 = uv1(:,in);
% ptrue = p(:,in);
% q2 = uv2(:,in);
% 
% q1(1:2,:) = q1(1:2,:);
% q2(1:2,:) = q2(1:2,:);
% 
% P1 = K*[R1 t1];
% P2 = K*[R2 t2];
% 
% precons = zeros(3,size(q1,2));
% 
%     %% here starts your code
% %    Q=[P1(3,:).*ptdraw(1,1)-P1(1,:);P1(3,:).*ptdraw(2,1)-P1(2,:);P2(3,:).*ptdraw(1,2)-P2(1,:);P2(3,:).*ptdraw(2,2)-P2(2,:)];
% %    [Uq Sq Vq]=svd(Q);
% %   X=Vq(end);
% %   for i=1:size(q1,2)
% %   precons(i,:)=;
% %   end
%        
% 
% %% visualization
% figure
% plot3(ptrue(1,:),ptrue(2,:),ptrue(3,:),'ro','MarkerSize',8);hold on;
% plot3(precons(1,:),precons(2,:),precons(3,:),'g+','MarkerSize',8);hold on;
% xlabel('x:(m)');
% ylabel('y:(m)');
% zlabel('z:(m)');
% title('Triangulation Simulation');
% legend({'Truth','Reconstruction'});
% set(gca,'FontName','Arial','FontSize',20);
% grid on;    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   END   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%