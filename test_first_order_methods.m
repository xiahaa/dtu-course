function test_first_order_methods
    clc;
    clear all;
    close all;
     
    skip_exe3 = 1;
    skip_solve = 1;
    skip_solve_trade_off = 1;
    
    if skip_exe3 == 0
        m = 5;
        n = 3;

        A = rand(m,n);
        b = rand(m,1);
        gamma = 0.4;
        AtA = A'*A;
        Atb = A'*b;

        %% PGM method
        [xPGM xPGMs] = solViaPGM(AtA,Atb,gamma,n);

        %% FPGM
        [xFPGM xFPGMs] = solViaFPGM(AtA, Atb, gamma,n);

        %% CVX solution
        cvx_begin
            variable x(n)
            minimize( 0.5*square_pos(norm(A*x-b,2)) + gamma*norm(x,1) )
        cvx_end

        disp(strcat('optimal solution from PGM:',mat2str(xPGM)));
        disp(strcat('optimal solution from FPGM:',mat2str(xFPGM)));
        disp(strcat('optimal solution from CVX:',mat2str(x)));
    end
    
    %% load l1regls.mat
    load('l1regls.mat');
    colnrm2 = sqrt(sum(A.^2,1));
    A = A ./ colnrm2;
    
    m = size(A,1);
    n = size(A,2);
    
    gamma = 0.4;
    AtA = A'*A;
    Atb = A'*b;
    
    if skip_solve == 0
        %% PGM method
        [xPGM xPGMs] = solViaPGM(AtA,Atb,gamma,n);

        %% FPGM
        [xFPGM xFPGMs] = solViaFPGM(AtA, Atb, gamma,n);

        %% CVX solution
        cvx_begin
            variable x(n)
            minimize( 0.5*square_pos(norm(A*x-b,2)) + gamma*norm(x,1) )
        cvx_end
        fstar = evaluateFunc(A,x,b,gamma);

        disp(strcat('optimal solution from PGM:',mat2str(xPGM)));
        disp(strcat('optimal solution from FPGM:',mat2str(xFPGM)));
        disp(strcat('optimal solution from CVX:',mat2str(x)));

        fPGMs = zeros(size(xPGMs,2),1);
        fFPGMs = zeros(size(xFPGMs,2),1);
%         for i=1:size(xPGMs,2)
        fPGMs = evaluateFunc(A,xPGMs,b,gamma);
%         end
%         for i=1:size(xFPGMs,2)
        fFPGMs = evaluateFunc(A,xFPGMs,b,gamma);
%         end
        save('solution.mat','xPGMs','xFPGMs','fPGMs','fFPGMs');
    else
        load('solution.mat');
    end
    %% drawing
    figure; 
    semilogy(fPGMs,'b');hold on; grid on;
    semilogy(fFPGMs,'r');
    legend('Proximal Gradient Method','Fast Proximal Gradient Method');
   
    %% trade-off curve
    if skip_solve_trade_off == 0
        maxGamma = max(abs(Atb));
        gammaSize = 200;
        gInterval = linspace(0,maxGamma,gammaSize);
        pairSols = zeros(gammaSize,2);
        for i = 1:1:numel(gInterval)
            gamma = gInterval(i);
            %% CVX solution
            cvx_begin
                variable x(n)
                minimize( 0.5*square_pos(norm(A*x-b,2)) + gamma*norm(x,1) )
            cvx_end
            pairSol(i,:) = [0.5*norm(A*x-b,2)^2, norm(x,1)];
        end
        save('solution_tradeoff.mat','pairSol','gInterval');
    else
        load('solution_tradeoff.mat');
    end
    
    figure;
    plot(pairSol(:,1),pairSol(:,2));
    title('Trade-off curve');
    xlabel('0.5*norm(Ax-b)^2');
    ylabel('one norm of x');
    grid on;
    
    figure;
    plot(gInterval',pairSol(:,1)+gInterval'.*pairSol(:,2));
    title('Cost');
    xlabel('gamma');
    ylabel('cost');
    grid on;
    
end

function [xsol xhistory]= solViaFPGM(AtA,Atb,gamma,n)
    L = norm(AtA).^2;
    t = 1/L;
    k = 0;
    TOLERANCE = 1e-6;
    MAXITER = 1e4;
    
    xk = zeros(n,1);
    xk_1 = xk;
    
    xhistory = [];
    
    while k <= MAXITER
        xhistory = [xhistory xk];
        if checkOptimalityConditions(AtA,Atb,xk,gamma,TOLERANCE) == 1
            break;
        end
        
        yk = xk + (k-1)/(k+2)*(xk-xk_1);
        %% gradient
        grad = computeGradient(AtA,Atb,yk);
        %% update 
        xnew = PGM(yk, t, grad, gamma);
    
        xk_1 = xk;
        xk = xnew;
        
        k = k+1;
    end
    xsol = xk;
end

function [xsol xhistory] = solViaPGM(AtA,Atb,gamma,n)
    L = norm(AtA).^2;
    xPGM = zeros(n,1);
    t = 1/L;
    k = 1;
    TOLERANCE = 1e-6;
    MAXITER = 1e4;
    xhistory = [];
    while k <= MAXITER
        xhistory = [xhistory xPGM];
        if checkOptimalityConditions(AtA,Atb,xPGM,gamma,TOLERANCE) == 1
            break;
        end
        %% gradient
        grad = computeGradient(AtA,Atb,xPGM);
        %% update 
        xnew = PGM(xPGM, t, grad, gamma);
        
        xPGM = xnew;
        
        k = k+1;
    end
    xsol = xPGM;
end

function flag = checkOptimalityConditions(AtA,Atb,x,gamma,tol)
    optimalityValue = AtA*x-Atb+gamma*sign(x);
    checkVal = norm(optimalityValue,1);
    if checkVal < tol
        flag = 1;
    else
        flag = 0;
    end
end

function fval = evaluateFunc(A,x,b,gamma)
    eval1 = A*x-repmat(b,1,size(x,2));
    colnrm2 = sum(eval1.^2,1);
    fval = 0.5 .* colnrm2 + gamma.*sum(abs(x),1);
%     fval = 0.5.*norm(A*x-b,2).^2+gamma*norm(x,1);
end

function grad = computeGradient(AtA,Atb,x)
    grad = AtA*x-Atb;
end

function xnew = PGM(x, t, grad, gamma)
    xeval = x-t*grad;
    h = gamma * t;
    
    idp = xeval >= h;
    idn = xeval <= -h;
    idz = abs(xeval) <= h;
    
    xnew = x;
    xnew(idp) = xeval(idp) - h;
    xnew(idn) = xeval(idn) + h;
    xnew(idz) = 0;    
end