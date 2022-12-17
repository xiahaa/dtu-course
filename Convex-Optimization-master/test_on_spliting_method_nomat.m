function test_on_spliting_method_nomat
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
%     D1 = spdiags([-ones(M-1,1),ones(M-1,1);0 1],[0,1],M,M);
%     D2 = spdiags([-ones(N-1,1),ones(N-1,1);0 1],[0,1],N,N);
%     D = vertcat(kron(speye(N),D1),kron(D2,speye(M)));
%     normA = powerIteration(gamma.*D);
    
    normA = 1;
    clear D D1 D2
    
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
    
    xk = zeros(M*N,1);
    yk = zeros(size(xk,1)*2,1);
    
    tv = @(x,gamma,M,N) gamma*sum(sqrt(sum(reshape(GetDx(x,M,N),[],2).^2)));
    tic
    maxIter = 10;
    k = 1;
    fs = [];
    while 1
        fval = evaluateCost(xk, b) + tv(xk,gamma,M,N);
        Dty = GetDty(yk, M, N);
        Dx = GetDx(xk,M,N);
        if checkOptimalityConditions(xk,b,gamma,Dty,Dx) == 1 || k > maxIter
            break;
        end
        grad = computeGradientFx(xk,b);
        xnew = proxg(xk - tau*(grad+gamma.*Dty),lb,ub);
        Dx = GetDx(2*xnew-xk,M,N);
        ynew = proxhstar(yk+sigma*gamma.*Dx,1);
        xk = rho*xnew + (1-rho)*xk;
        yk = rho*ynew + (1-rho)*yk;
        k=k+1;
        fs = [fs;fval];
    end
    toc
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

function Dx = GetDx(xk,M,N)
    Dx = zeros(M*N*2,1);
    for i = 2:M-1
        for j = 2:N-1
%             if i < M
            Dx((j-1)*M+i) = xk((j-1)*M+i+1) - xk((j-1)*M+i);
%             end
%             if j < N
            Dx((j-1)*M+i+M*N) = xk((j)*M+i) - xk((j-1)*M+i);
%             end
        end
    end
end        

function Dty = GetDty(y, M, N)
    % D is 2NM x MN, so Dt is MN x 2MN, y is 2MN x 1
    Dty = zeros(M*N,1);
    
    for i = 2:M-1
        for j = 2:N-1
            % up
%             if (i - 1) > 0
            Dty((j-1)*M+i) = Dty((j-1)*M+i)+y((j-1)*M+i-1);
%             end
            % down
%             if (i + 1) <= M
            Dty((j-1)*M+i) = Dty((j-1)*M+i)-y((j-1)*M+i);
%             end
            % left
%             if (j - 1) > 0
            Dty((j-1)*M+i) = Dty((j-1)*M+i)+y((j-2)*M+M*N+i);
%             end
            % right
%             if (j + 1) <= N
            Dty((j-1)*M+i) = Dty((j-1)*M+i)-y((j-1)*M+M*N+i);
%             end
        end
    end
end

function flag = checkOptimalityConditions(x,b,gamma,Dty,Dx)
    cond1 = norm(-computeGradientFx(x,b)-gamma.*Dty,1);
    cond2 = norm(-Dx,1);
    
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