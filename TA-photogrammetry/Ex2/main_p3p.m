clc;close all;clear all;
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

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
xlabel('x:(m)','FontName','Arial','FontSize',15);
ylabel('y:(m)','FontName','Arial','FontSize',15);
zlabel('z:(m)','FontName','Arial','FontSize',15);
grid on;
title('P3P Simulation','FontName','Arial','FontSize',15);

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





