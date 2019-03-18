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


% I = drawlines(I,q,[[1 2];[3 4];[1 4];[2 3]]);

function ph = tohomogeneous(p)
    ph = [p;ones(1,size(p,2))];
end

function p = fromhomogeneous(ph)
    p = ph(1:end-1,:) ./ ph(end,:);
end

function rot = Rxyz(roll, pitch, yaw)
% function rot = Rxyz(roll, pitch, yaw)
%
% Create a rotation matrix corresponding to the given Euler angles. The
% rotation order is assumed to be fistly along z-axis, then y-axis, and
% finally x-axis.
%   Inputs:
%       roll: rotation angle along x-axis, rad.
%       pitch: rotation angle along y-axis, rad.
%       yaw: rotation angle along z-axis, rad.
%   Outputs:
%       rot: Rotation Matrix
%
% Author: xiahaa@space.dtu.dk
% Disclaimer: This code comes with no guarantee at all and its author
%   is not liable for any damage that its utilization may cause.
    rot1 = [cos(yaw) sin(yaw) 0;-sin(yaw) cos(yaw) 0; 0 0 1];
    rot2 = [cos(pitch) 0 -sin(pitch);0 1 0; sin(pitch) 0 cos(pitch)];
    rot3 = [1 0 0;0 cos(roll) sin(roll);0 -sin(roll) cos(roll)];

    rot = rot3*rot2*rot1;
    rot = rot';
end

