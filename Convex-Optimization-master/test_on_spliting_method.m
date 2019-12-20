function test_on_spliting_method
    clear all;
    close all;
    clc;
    
    % Get cameraman test image
    fname = fullfile(matlabroot,'toolbox','images',...
                     'imdata','cameraman.tif');
    X = double(imread(fname))/255;
    [M,N] = size(X);
    % Generate pseudo-random noisy image
    b = X(:) + 0.05*randn(M*N,1);

    m = M;n = N;
    gamma = 0.1;
    lb = 0; ub = 1;
    
    %% data preparation
    % Finite-difference matrix with Neumann boundary conditions
    D1 = spdiags([-ones(M-1,1),ones(M-1,1);0 1],[0,1],M,M);
    D2 = spdiags([-ones(N-1,1),ones(N-1,1);0 1],[0,1],N,N);
    D = vertcat(kron(speye(N),D1),kron(D2,speye(M)));
    normA = powerIteration(gamma.*D);
    
    L = 1;
    sigma = 0.8;
    tau = 1/(L*0.5+sigma*normA) * 0.9;
%     tau = 0.1;
%     L = 1;
%     sigma = ((1/tau-L/2)/(normA))*0.5;
    
    delta = 2-L/2/(1/tau-sigma*normA);
    if delta < 1
        delta = 1;
    elseif delta > 2
        delta = 2;
    end
    rho = 1;
    
    xk = zeros(numel(b),1);
    yk = zeros(size(D,1),1);
    
    tv = @(x,gamma) gamma*sum(sqrt(sum(reshape(D*x,[],2).^2)));
    
    maxIter = 1e2;
    k = 1;
    fs = [];
    while 1
        fval = evaluateCost(xk, b) + tv(xk,gamma);
        
        if checkOptimalityConditions(xk,yk,b,D,gamma) == 1 || k > maxIter
            break;
        end
        grad = computeGradientFx(xk,b);
        xnew = proxg(xk - tau*(grad+gamma.*D'*yk),lb,ub);
        ynew = proxhstar(yk+sigma*gamma.*D*(2*xnew-xk),1);
        xk = rho*xnew + (1-rho)*xk;
        yk = rho*ynew + (1-rho)*yk;
        k=k+1;
        fs = [fs;fval];
    end
    noisyImage = reshape(b,M,N);
    recI = reshape(xk,M,N);
    imcombine = double(zeros(M,N*3));
    imcombine(1:M,1:N) = X;
    imcombine(1:M,N+1:N*2) = noisyImage;
    imcombine(1:M,2*N+1:N*3) = recI;
    imshow(imcombine,[]);
    figure
    plot(fs);
end

function feval = powerIteration(A)
    maxIter = 1000;
    AtA = A'*A;
    r = rand(size(AtA,1),1);
    for i=1:maxIter
        r1 = AtA*r;
        r = r1./norm(r1);
    end
    feval = r'*AtA*r;
end

function feval = fun_A(x,mode)
    %% assume A = I;
    if mode == 1
        feval = x;
    elseif mode == 2
        feval = x;
    end
end

function flag = checkOptimalityConditions(x,y,b,D,gamma)
    cond1 = norm(-computeGradientFx(x,b)-gamma.*D'*y,1);
    cond2 = norm(-D*x,1);
    
    tol = 1e-4;
    
    if cond1 < tol && cond2 < tol
        flag = 1;
    else
        flag = 0;
    end
end

function fx = evaluateCost(x, b)
    Ax = fun_A(x,1);
    fx = 0.5*norm(Ax-b,2).^2;
end

function grad = computeGradientFx(x,b)
    AtAx = fun_A(fun_A(x,1),2);
    Atb = fun_A(b,2);
    grad = AtAx-Atb;
end

function proxGrad = proxg(x,lb,ub)
    idp = x >= ub;
    idn = x <= lb;
    proxGrad = x;
    proxGrad(idp) = ub;
    proxGrad(idn) = lb;
end

function proxGrad = proxhstar(x,ub)
    proxGrad = x;
%     if norm(x)>ub
%         proxGrad = x./norm(x).*ub;
%     end
    id = abs(x)>ub;
    proxGrad(id)=sign(x(id))*ub;
end