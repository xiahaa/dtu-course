close all;
u = [[0 0];[2 2];[4 3]];
Omega(1,:,:) = [1 0.5;0.5 1];
Omega(2,:,:) = [1 -0.7;-0.7 1];
Omega(3,:,:) = [0.2 0.1;0.1 0.5];
w = [0.5;0.2;0.3];
ii = 1;
probs = [];
for x1 = -4:0.1:6
    jj = 1;
    for x2 = -4:0.1:6
        probs(ii,jj) = multigaussian([x1 x2]', u, Omega,w);
        jj = jj + 1;
    end
    ii = ii + 1;
end

figure
surf([-4:0.1:6]',[-4:0.1:6]',probs);
