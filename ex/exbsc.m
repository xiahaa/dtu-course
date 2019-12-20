clc; clear; close all;

ex4bis

function ex1 
    n = 200;
    h = 1/n;
    D= 0.1;
    aux = zeros(n-1,1);
    aux(1) = -2;
    aux(2) = 1;
    L = (toeplitz(aux));
    fun = @(t,u) (D/(h^2).*(L*u));
    
    x = 1/n*(1:n-1)';
    u0 = zeros(n-1,1);
    u0(x>1/3&x<2/3) = 1;
    
    figure(1);
    plot([0;x;1],[0;u0;0],'r-','LineWidth',3);
    axis([0,1,0,1.1]);
    set(gca,'FontSize',20);
    
    tspan = [0:0.02:1];
    [time,u] = ode15s(fun,tspan,u0);
    
    figure(2);
    for j = 1:length(time)
        plot([0;x;1],[0;u(j,:)';0],'k-','LineWidth',3);
        text(0.2,2,['t=' num2str(time(j))],'FontSize',15);
        axis([0,1,0,2.2]);
        pause(0.01);
    end
    
    % backward
    fun2 = @(t,u) (-D/(h^2).*(L*u));
    ufinal = u(end,:)';
    ufinal = ufinal + 0.01*rand(size(ufinal));
    [time2,u2]=ode15s(fun2,tspan,ufinal);
    
    figure(3);
    for j = 1:length(time2)
        plot([0;x;1],[0;u2(j,:)';0],'k-','LineWidth',3);
        text(0.2,2,['t=' num2str(time2(j))],'FontSize',15);
        axis([0,1,0,2.2]);
        pause(1);
    end
    
end


function ex2
    n = 200;
   
    x0 = randn(n,1);
    [L1,L2]=estab_l1l2(n);
    sigma = 1;
    figure
    m = 10;
    colormap = 'jet';
    for i = 1:m
        w = sigma.*randn(n,1);
        x = x0 + L1\w;
        plot(x);hold on;
    end
    
    figure
    m = 10;
    for i = 1:m
        w = sigma.*randn(n,1);
        x = x0 + L2\w;
        plot(x);hold on;
    end
    
    [~,L2]=estab_l1l2(n-1);
    L2 = (n-1)^2.*L2;
    sigma = 1;
    figure
    m = 10;
    for i = 1:m
        w = sigma.*randn(n-1,1);
        x = L2\[w];
        plot([0;x;0]);hold on;
    end
    
    %% task 3-4
    n = 200;
    [~,L2]=estab_l1l2(n-1);
    L2 = (n-1)^2*L2;
    N = (n-1);
    L1 = kron(speye(N),L2);
    L2 = kron(L2,speye(N));
%     spy(L1);
%     spy(L2);
    L = L1+L2;
    lambda = 0.05;
    for i = 1:2
        lambda=lambda * 10^(i-1);
        M = L-1/lambda^2.*speye(N*N);
        W = randn(N*N,1);
        X = M\W;
        X = reshape(X,n-1,n-1);
        X_0 = zeros(n+1,n+1);
        X_0(2:end-1,2:end-1)=X;
        figure
        imagesc(X_0);
        title(['$\lambda$=',num2str(lambda)],'Interpreter','latex');
    end
    
    %% adding structure
    n = 100;
    [L1,L2]=estab_l1l2(n-1);
    L1 = (n-1)*L1;
    L2 = (n-1)^2.*L2;
    sigma = 1;
    figure
    m = 10;
    id = randperm(n-1,3);
    for i = 1:m
        w = sigma.*randn(n-1,1);
        w(id) = w(id).*10;
        x = L1\[w];
        plot([0;x;0]);hold on;
        plot(id+1,x(id),'ko');
    end
    
    figure
    for i = 1:m
        w = sigma.*randn(n-1,1);
        w(id) = w(id).*10;
        x = L2\[w];
        plot([0;x;0]);hold on;
        plot(id+1,x(id),'ko');
    end
    
    %% hypermodel
    n = 100;
    m = 10;
    [L1,L2]=estab_l1l2(n-1);
    L1 = (n-1)*L1;
    L2 = (n-1)^2.*L2;
    theta = 10*randn(n-1,1)';
    Dtheta = diag(theta);
    figure
    subplot(1,2,1);
    for i = 1:m
        w = Dtheta*randn(n-1,1);
        x2 = L2\w;
        x1 = L1\w;
        subplot(1,2,1); 
        plot([0;x1;0]);hold on;
        subplot(1,2,2); 
        plot([0;x2;0]);hold on;
    end
end

function ex3
    n = 100; h = 1/n; N = n-1;
    [L1,L2] = estab_l1l2(N);
    %% laplacian
    L1 = 1/h .* L1;
    L2 = 1/h^2 .* L2;
    
    %%
    bk = [25;35;50;80];
    A = zeros(length(bk),N);
    for i = 1:length(bk)
        A(i,bk(i)) = 1;
    end
    b = [0.5,1,0.2,2]';
    sigma = 1;
    
    figure
    colormap = 'jet';
    lambda = 1;
    for i = 1:10
        gamma = 0.01 * 10^(i-1);
        %% one-dimensional Whittle-Matern prior
        M = L2-1/lambda^2.*speye(N);
        %% ls
        lhs = [1/gamma*M;1/sigma*A];
        rhs = [zeros(N,1);1/sigma*b];
        x = lhs \ rhs;
        plot([0;x;0],'LineWidth',1.5);hold on;
    end
    
end

function ex4
    dfun = @(t) (12/sqrt(pi).*exp(-(6.*t-3).^2));
    fun = @(t) (1+erf(6.*t-3));
    
    n = 100;
    h = 1/n;
    t = 0:h:1;
    
    sigma = 0.01;
    gamma = 500;
    
    y = fun(t);
    y = y + sigma.*randn(size(y));
    dy = dfun(t);
    dy1 = [diff(y)./h,0];
    
    %% task 1
    figure
    plot(y,'LineWidth',1.5);hold on; grid on;
    plot(dy,'LineWidth',1.5);
    plot(dy1,'LineWidth',1.5);
    legend({'y','dy-analytical','dy-numerical'});
    
    %% task 2
    A = toeplitz(ones(n-1,1),[1,zeros(1,n-2)]);
    A = A.*h;
    
    [L1,L2] = estab_l1l2(n-1);
    L1 = L1/h;
    L2 = L2/h^2;
    
    lambda = 0.2;
    M = L2 - 1/lambda^2.*speye(n-1);
    b = y(2:end-1)';
    lhs = [1/sigma*A;1/gamma*M];
    rhs = [1/sigma*b;zeros(n-1,1)];
    
    %% this is mean
    x = lhs \ rhs;
    
    %% covariance.
    D = gamma^2.*speye(n-1);
    Sigma = sigma^2.*speye(n-1);
    C = D-D*A'*inv(A*D*A'+Sigma)*A*D;
    cij = sqrt(diag(C));
    
    x = [0;x;0];
    cij = [0;cij;0];
    
    xlow2 = x -2*cij;
    xhigh2 = x+2*cij;
    xlow1 = x - cij;
    xhigh1 = x+ cij;
    
    figure
    fill([t';t(n+1:-1:1)'],[xlow2;xhigh2(n+1:-1:1)],[0.95,0.95,1])
    hold on
    fill([t';t(n+1:-1:1)'],[xlow1;xhigh1(n+1:-1:1)],[1,0.95,0.95])
    set(gca,'FontSize',20)
    
    plot(t,dy,'k-','LineWidth',2);hold on;
    plot(t,x,'r-','LineWidth',2);
    

    
%     plot(t,x-2*cij,'k--','LineWidth',3);
%     plot(t,x+2*cij,'k--','LineWidth',3);
end

function test_cg_sol
    A = [1 2;3 4];
    x = [1.2;3.4];
    b = A*x;
    
    fun=@(x) (A*x);
    res=@(bhat) (b-bhat);
    
    xsol = cg_sol(fun, res, A);
    disp(xsol-x)
end

function y = fun_a(t, s)
    % the width parameter
    w = 0.03; 
    % inverse of the width parameter.
    k = 1/w;
    fun_j = @(x) (besselj(1,x));
    y = zeros(1,length(s));
    id1 = s ~= t;
    id2 = s == t;
    y(id2) = 1;
    y(id1) = (fun_j(k.*abs(t-s))./(k.*abs(t-s))).^2;
end

function An = est_An(N, tj)
    k = 1:N;
    s = (k - 0.5)/N;
    An = zeros(length(tj),N);
    for j = 1:length(tj)
        An(j,:) = fun_a(tj(j), s)/N;
%         plot(An(j,:));
    end
end

function ex4bis
    %% define relevant function handler
    m = 50;
    sigma = 0.0003;
    j = 1:m;
    tj = (j-0.5)/m;
    
    %% linspace 101, h = 0.1
    N = 101;
    An = est_An(N, tj);
    k = 1:N;
    id = k > N/2 & k < N*2/3;
    xn = zeros(N,1); xn(id) = 1;
    
    b = An * xn;
    bn = b + sigma*randn(size(b));
    
    figure
    plot(tj, b, 'b-', 'LineWidth',2);hold on; grid on;
    plot(tj, bn, 'r-', 'LineWidth',2);
    legend({'true','noisy'});
    set(gca,'FontSize', 20);
    
    n = 200;
    A = est_An(n, tj);
    
    [L1,~] = estab_l1l2(n);
    L = L1*n;
    Linv = inv(L);
    %% posterior:
    % solving using IAS algorithm
    theta_star = 1e4*ones(n,1);% prior
    beta = 3/2+1e-2;
    eta = beta - 3/2;
    theta = theta_star;
    old_w = [];
    while true
        %% known theta, update w
        D = diag([1./sqrt(theta)]);
        lhs = [1/sigma*A*Linv;D];
        rhs = [1/sigma*bn;zeros(n,1)];
        w = lhs \ rhs;
        %% know w, update theta
        theta_k = theta_star / 2 .* (eta + sqrt(eta^2+2*w.^2./theta_star));
        %% exist
        if norm(theta_k - theta)<1e-3 && (~isempty(old_w) || norm(old_w-w) < 1e-3)
            break;
        end
        disp('iter')
        theta = theta_k;
        old_w = w;
    end
    x = Linv*w;
    figure
    plot(linspace(0,1,length(xn)),xn,'b-','LineWidth',2);hold on;grid on;
    plot(linspace(0,1,length(x)),x,'r-','LineWidth',2);
    legend({'real','inverse'});
    set(gca,'FontSize',20);
end

function x = cg_sol(f, r, A)
   x0 = zeros(2,1);
   r0 = r(f(x0));
   p0 = r0;
   
   while true
       %
       a = norm(r0)^2/(p0'*A*p0);
       x1 = x0 + a*p0;
       r1 = r0 - a*A*p0;
       b = norm(r1)^2/norm(r0)^2;
       p1 = r1+b*p0;
       
       if norm(x1-x0)<1e-10
           break;
       end
       
       x0 = x1;
       r0 = r1;
       p0 = p1;
   end
   x = x0;
end
