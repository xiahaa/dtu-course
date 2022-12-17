clc;clear all;close all;
M = 3;
N = 4;

% Finite-difference matrix with Neumann boundary conditions
D1 = spdiags([-ones(M-1,1),ones(M-1,1);0 1],[0,1],M,M);
D2 = spdiags([-ones(N-1,1),ones(N-1,1);0 1],[0,1],N,N);
D = vertcat(kron(speye(N),D1),kron(D2,speye(M)));

x = randn(M*N,1);
gamma = 1;
tv1 = @(x,gamma) gamma*sum(sqrt(sum(reshape(D*x,[],2).^2)));

tv2 = @(x,gamma,M,N) gamma*sum(sqrt(sum(reshape(GetDx(x,M,N),[],2).^2)));

tv1(x,gamma)
tv2(x,gamma,M,N)

y = D*x;

D'*y

GetDty(y, M, N)


function Dx = GetDx(xk,M,N)
    Dx = zeros(M*N*2,1);
    for i = 1:M
        for j = 1:N
            if i < M
                Dx((j-1)*M+i) = xk((j-1)*M+i+1) - xk((j-1)*M+i);
            end
            if j < N
                Dx((j-1)*M+i+M*N) = xk((j)*M+i) - xk((j-1)*M+i);
            end
        end
    end
end  

function Dty = GetDty(y, M, N)
    % D is 2NM x MN, so Dt is MN x 2MN, y is 2MN x 1
    Dty = zeros(M*N,1);
    
    for i = 1:M
        for j = 1:N
            % up
            if (i - 1) > 0
                Dty((j-1)*M+i) = Dty((j-1)*M+i)+y((j-1)*M+i-1);
            end
            % down
            if (i + 1) <= M
                Dty((j-1)*M+i) = Dty((j-1)*M+i)-y((j-1)*M+i);
            end
            % left
            if (j - 1) > 0
                Dty((j-1)*M+i) = Dty((j-1)*M+i)+y((j-2)*M+M*N+i);
            end
            % right
            if (j + 1) <= N
                Dty((j-1)*M+i) = Dty((j-1)*M+i)-y((j-1)*M+M*N+i);
            end
        end
    end
end