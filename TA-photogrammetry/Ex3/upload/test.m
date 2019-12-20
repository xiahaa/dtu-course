baseDir = './images/';% use you own directory
buildingScene = imageDatastore(baseDir);
numImages = numel(buildingScene.Files);

load(strcat(baseDir,'viff.xy'));
x = viff(1:end,1:2:72)'; % pull out x coord of all tracks
y = viff(1:end,2:2:72)'; % pull out y coord of all tracks

% visualization
num = 0;
for n = 1:numImages-1
    im1 = readimage(buildingScene, n);
    imshow(im1); hold on;
    id = x(n,:) ~= -1 & y(n,:) ~= -1;
    plot(x(n,id),y(n,id),'go');
    num = num + sum(id);
    hold off;
    pause(0.1);
end

% load projection matrices
load(strcat(baseDir,'dino_Ps.mat'));

ptcloud = zeros(3,num);
k = 1;
for i = 1:size(x,1)-1
% tracked features
    id = x(i,:) ~= -1 & y(i,:) ~= -1 & x(i+1,:) ~= -1 & y(i+1,:) ~= -1;
    q1 = [x(i,id);y(i,id);];
    q2 = [x(i+1,id);y(i+1,id);];
    P1 = P{i};
    P2 = P{i+1};
    
    precons = zeros(3,size(q1,2));
    for j = 1:size(q1,2)
        % your code starts
        precons(:,j) = triangulationMidpoint(q1(:,j),P1,q2(:,j),P2);
    end

    ptcloud(:,k:k+size(q1,2)-1) = precons;
    k = k + size(q1,2);
end

figure
plot3(ptcloud(1,:),ptcloud(2,:),ptcloud(3,:),'k.','MarkerSize',10);hold on;
grid on;
axis equal;
view(3);
for i = 1:size(x,1)
    P1 = Pfig;
    [K, R, t, c] = decomposeP(P1);
    plotCamera('Location',c,'Orientation',R,'Opacity',0,'Color',[0 1 0],'Size',0.05);
end