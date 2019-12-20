function test_on_gdop
    clear all;
    clc;
    N = 10000;
    dop1 = zeros(N,1);
    dop2 = zeros(N,1);
    t1 = zeros(N,1);
    t2 = zeros(N,1);
    
    for j = 1:N
        j
        m = randi(50,1) + 4;
%         m = 10;
        H = randn(m,3);
        Hn = sqrt(H(:,1).^2+H(:,2).^2+H(:,3).^2);
        for i = 1:size(H,1)
            H(i,1:3) = H(i,1:3) ./ Hn(i);
        end
        Ha = [H ones(size(H,1),1)];
        M = Ha'*Ha;
    
        %% 1
        tic
        dop = dopsss(M);
        t1(j)=toc;
        dop1(j) = dop;
    
        %% 2
        tic
        dop = dopschur(M);
        t2(j)=toc;
        dop2(j) = dop;
    end
    figure(1);
    plot(dop1,'r');hold on;
    plot(dop2,'b');grid on;
    
%     figure(2);
%     plot(t1,'r');hold on;
%     plot(t2,'b');grid on;
    disp(mean(t1));
    disp(mean(t2));
end

function dop = dopsss(M)
    Minv = inv(M);
    dop = trace(Minv(1:3,1:3));
end

function dop = dopschur(M)
    A = M(1:3,1:3);B = M(1:3,4);
    S = A - B*B'/M(4,4);
    e = eig(S);
    dop = (1/e(1)+1/e(2)+1/e(3));
end

function dop = dopeee(M)
    e1 = eig(M);
    A = M(1:3,1:3);h = M(1:3,4);
    [evec,eval] = eig(A);
    hh = evec' * h;
    e = diag(eval);
    hhh = hh(1)*hh(1)*1/e(1) + hh(2)*hh(2)*1/e(2) + hh(3)*hh(3)*1/e(3);
    z = 1 / (M(4,4)-hhh);
    dop = sum(1./e1)-z;
end

