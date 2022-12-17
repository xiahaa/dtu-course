clc;clear all;close all;

% Get cameraman test image
fname = fullfile(matlabroot,'toolbox','images',...
                 'imdata','cameraman.tif');
X = double(imread(fname))/255;
[M,N] = size(X);
% Generate pseudo-random noisy image
b = X(:) + 0.05*randn(M*N,1);

m = M;n = N;

u = zeros(M*N,1);
sigma_x = zeros(M*N,1);
sigma_y = zeros(M*N,1);

miu = 1;
beta = 1;

D1 = spdiags([-ones(M-1,1),ones(M-1,1);0 1],[0,1],M,M);
D2 = spdiags([-ones(N-1,1),ones(N-1,1);0 1],[0,1],N,N);
D1 = kron(speye(N),D1);
D2 = kron(D2,speye(M));
D = vertcat(D1,D2);

wx = shrink(D1*u+sigma_x, 1/beta);
wy = shrink(D2*u+sigma_y, 1/beta);
eta = zeros(M*N,1);

maxIter = 100;
A = speye(M*N);

tv = @(x,gamma) gamma*sum(sqrt(sum(reshape(D*x,[],2).^2)));

tmp = beta * D1'*D1 + beta * D2'*D2 + miu * A'*A;
uold = zeros(size(u));
costs = [];
tic
for i = 1:maxIter
    costs(i) = 0.5*norm(A*u-b,2).^2 + tv(u,beta);
    d = D1'*sigma_x + D2'*sigma_y + A'*eta + beta * D1' * (D1 * u - wx) + beta * D2' * (D2 * u - wy) + miu * A' * (A*u-b);
    alpha_k = d'*d/(d'*tmp*d);
    uold = u;
    u = u - alpha_k * d;
    wx = shrink(D1*u+sigma_x/beta, 1/beta);
    wy = shrink(D2*u+sigma_y/beta, 1/beta);
    
    ee = norm(u-uold,'fro') / norm(u,'fro');
    if (ee < 1e-4)
        if (ee < 1e-5)
            break;
        end
        sigam_x = sigma_x + beta * (D1*u-wx);
        sigam_y = sigma_y + beta * (D2*u-wy);
        eta = eta + miu * (A*u-b);
        beta = beta * 2;
        if beta > 2^9
            beta = 2^9;
        end
        miu = 2 * miu;
        if miu > 2^13
            miu = 2^13;
        end
    end
end
toc

i

noisyImage = reshape(b,M,N);
recI = reshape(u,M,N);
imcombine = double(zeros(M,N*3));
imcombine(1:M,1:N) = X;
imcombine(1:M,N+1:N*2) = noisyImage;
imcombine(1:M,2*N+1:N*3) = recI;
imshow(imcombine,[]);
figure
plot(costs);

function res = shrink(u, alpha)
    res = zeros(numel(u),1);
    ind = abs(u) > alpha;
    res(ind) = sign(u(ind)).*(abs(u(ind))-alpha);
end


function Dx = GetDx(u,M,N)
    Dx = zeros(M*N,1);
    for i = 1:M
        for j = 1:N
            if j < N
                Dx((j-1)*M+i) = u((j)*M+i) - u((j-1)*M+i);
            end
        end
    end
end

function Dy = GetDy(u,M,N)
    Dy = zeros(M*N,1);
    for i = 1:M
        for j = 1:N
            if i < M
                Dy((j-1)*M+i) = u((j-1)*M+i+1) - u((j-1)*M+i);
            end
        end
    end
end