function test_fast_dop_computation
    n = 20;
    dim = 3;
    random_vecs = gen_ran_vecs(n, dim);
    p = v_normalize(random_vecs);
    
    runcnt = 100000;
    Hs = cell(runcnt,1);
    for k = 1:runcnt
        m = randi([4,20],1);
        sel = randperm(n,m);
        A = p(sel',:);
        H = [A ones(numel(sel),1)];
        Hs{k} = H;
    end
    
    t1s = zeros(runcnt,1);
    t2s = zeros(runcnt,1);
    t3s = zeros(runcnt,1);
    t4s = zeros(runcnt,1);
    t5s = zeros(runcnt,1);

    
    dop1 = zeros(runcnt,1);
    dop2 = zeros(runcnt,1);
    dop3 = zeros(runcnt,1);
    dop4 = zeros(runcnt,1);
    dop5 = zeros(runcnt,1);

    for k = 1:runcnt
        H = Hs{k};
        tic;
        dop1(k) = dop_via_inverse(H);
        t1s(k) = toc;
        
        %% 2
        tic;
        dop2(k) = dop_via_cl(H);
        t2s(k) = toc;
        
        %% 3
        tic;
        dop3(k) = eig_decomp(H);
        t3s(k) = toc;
        
        %% 4
        tic;
        dop4(k) = trial(H);
        t4s(k) = toc;
 
        %% 5
        tic;
        dop5(k) = fast_eig_decomp(H);
        t5s(k) = toc;
    end
    
    
    
    mean(t1s)
    mean(t2s)
    mean(t3s)
    mean(t4s)
    mean(t5s)
    
    sum(t1s)
    sum(t2s)
    sum(t3s)
    sum(t4s)
    sum(t5s)
    
    plot(dop1,'r');hold on; grid on;
    plot(dop2,'b--');
    plot(dop3,'m-o');
    
end

function dop = eig_decomp2(H)
%     S = svd(H);
%     E = S.^2;
%     dop = sqrt(1/E(1) + 1/E(2) + 1/E(3) + 1/E(4));
    M = H'*H;
%     tic
    a = M(1,1);b = M(1,2);c = M(1,3);d = M(1,4);
    e = M(2,2);f = M(2,3);g = M(2,4);
    h = M(3,3);i = M(3,4);
    j = M(4,4);

    bb = b*b;
    cc = c*c;
    dd = d*d;
    ff = f*f;
    gg = g*g;
    ii = i*i;
    ae = a*e;
    eh = e*h;
    ah = a*h;
    aj = a*j;
    ej = e*j;
    hj = h*j;

    fgi = f*g*i;
    bcf = b*c*f;
    bdg = b*d*g;
    ehj = e*h*j;
    
%     p1 = 1;%%lambda^4 
%     p2 = -e - h - j - a;%%- e*lambda^3 - h*lambda^3 - j*lambda^3 - a*lambda^3 
%     p3 = -bb - cc - dd - ff - gg - ii + ae + ah + aj + eh + ej + hj;
    p4 = a*ff+ a*gg+ cc*e + dd*e+ a*ii+ bb*h + bb*j+ dd*h + cc*j+ e*ii + ...
        gg*h+ ff*j - 2*fgi- ehj- 2*bcf- 2*bdg- ae*h- ae*j- 2*c*d*i- ah*j;

    p5 = cc*gg + dd*ff + bb*ii - ae*ii - ah*gg - aj*ff - dd*eh - cc*ej - ...
        2*c*d*f*g + 2*bcf*j - 2*b*c*g*i - 2*b*d*f*i + 2*bdg*h + 2*c*d*e*i ...
        + 2*a*fgi - bb*hj + a*ehj;
  
    dop = sqrt(-p4/p5);
%     toc
    

end

function [p3,p4] = coeff3(M)
    M12_M21=M(1,2)*M(2,1);
    M11_M22=M(1,1)*M(2,2);
    M13_M31=M(1,3)*M(3,1);
    M23_M32=M(2,3)*M(3,2);
    
    % M1_2*M2_1 - M1_1*M2_2 - M1_1*M3_3 + M1_3*M3_1 - M2_2*M3_3 + M2_3*M3_2, 
    p3 = M12_M21 - M11_M22 - M(1,1)*M(3,3) + M13_M31 - M(2,2)*M(3,3) + M23_M32;
    
    % M1_1*M2_2*M3_3 - M1_1*M2_3*M3_2 - M1_2*M2_1*M3_3 + M1_2*M2_3*M3_1 + M1_3*M2_1*M3_2 - M1_3*M2_2*M3_1, 
    p4 = M11_M22*M(3,3)-M(1,1)*M23_M32-M12_M21*M(3,3)+M(1,2)*M(2,3)*M(3,1)+M(1,3)*M(2,1)*M(3,2)-M13_M31*M(2,2);
end

function p = coeff4(M)
    n = size(M,1);
    p = zeros(1,n+1);
    p(1) = 1;
    C = M;
    p(2) = -trace(C);
    for k = 2:n
        C(1,1)=C(1,1)+p(k);
        C(2,2)=C(2,2)+p(k);
        C(3,3)=C(3,3)+p(k);
        C = M*C;
        p(k+1)=-trace(C)/k;
    end
end

function dop = trial(H)
    M = H'*H;
    p1 = coeff4(M);
    S1 = M(1:3,1:3) - M(1:3,4)*M(4,1:3)./M(4,4);
    p2 = coeff4(S1);
    dop = sqrt(-p1(end-1)/p1(end));
    pdop = sqrt(-p2(end-1)/p2(end));
    tdop = sqrt(-p1(end-1)/p1(end)+p2(end-1)/p2(end));
end

function dop = fast_eig_decomp(H)
    M = H'*H;

    %% gdop
    a = M(1,1);b = M(1,2);c = M(1,3);d = M(1,4);
    e = M(2,2);f = M(2,3);g = M(2,4);
    h = M(3,3);i = M(3,4);
    j = M(4,4);

    bb = b*b;
    cc = c*c;
    dd = d*d;
    ff = f*f;
    gg = g*g;
    ii = i*i;
    ae = a*e;
    eh = e*h;
    ah = a*h;
    aj = a*j;
    ej = e*j;
    hj = h*j;

    fgi = f*g*i;
    bcf = b*c*f;
    bdg = b*d*g;
    ehj = e*h*j;
    
    p4 = a*ff+ a*gg+ cc*e + dd*e+ a*ii+ bb*h + bb*j+ dd*h + cc*j+ e*ii + ...
        gg*h+ ff*j - 2*fgi- ehj- 2*bcf- 2*bdg- ae*h- ae*j- 2*c*d*i- ah*j;

    p5 = cc*gg + dd*ff + bb*ii - ae*ii - ah*gg - aj*ff - dd*eh - cc*ej - ...
        2*c*d*f*g + 2*bcf*j - 2*b*c*g*i - 2*b*d*f*i + 2*bdg*h + 2*c*d*e*i ...
        + 2*a*fgi - bb*hj + a*ehj;
    dop = sqrt(-p4/p5);
%     p = coeff4(M);
%     A = H(:,1:3)'*H(:,1:3);
%     s1 = sum(H(:,1));s2 = sum(H(:,2));s3 = sum(H(:,3));
%     s1s1 = s1*s1;
%     s1s2 = s1*s2;
%     s1s3 = s1*s3;
%     s2s3 = s2*s3;
%     s2s2 = s2*s2;
%     s3s3 = s3*s3;
%     N = size(H,1);
%     tic

%     [e1,v1] = eig(M(1:3,1:3));
    S1 = M(1:3,1:3) - M(1:3,4)*M(4,1:3)./M(4,4);% M(1:3,1:3) - [s1s1 s1s2 s1s3;s1s2 s2s2 s2s3;s1s3 s2s3 s3s3]./N;%
%     tic
%     p = coeff4(S1);
%     ps3 = p(end-1);ps4 = p(end);
%     toc
%     tic
    [ps3,ps4] = coeff3(S1);
%     toc
    pdop = sqrt(-ps3/ps4);
    tdop = sqrt(-p4/p5+ps3/ps4);
%     S2 = M(4,4) - M(4,1:3)*e1*diag([1/v1(1,1),1/v1(2,2),1/v1(3,3)])*e1'*M(1:3,4);
%     e2 = eig(S1);     
%     pdop = sqrt(1/e2(1)+1/e2(2)+1/e2(3));
%     tdop = sqrt(dop*dop-pdop*pdop);

%     dop = sqrt(1/e2(1)+1/e2(2)+1/e2(3) + 1/S2);
%     toc
%     tic
%     Minv = inv(M);
%     toc
end

function dop = eig_decomp(H)
    E = eig(H'*H);
    dop = sqrt(1/E(1) + 1/E(2) + 1/E(3) + 1/E(4));
end

function dop = dop_via_cl(H)
    M = H'*H;
    h1 = trace(M);
    h2 = trace(M*M);
    h3 = trace(M*M*M);
    h4 = det(M);
    dop = sqrt((0.5*h1*h1*h1-1.5*h1*h2+h3)/(3*h4));
end

function dop = dop_via_inverse(H)
    Minv = inv(H'*H);
    dop = sqrt(trace(Minv));
    pdop = sqrt(trace(Minv(1:3,1:3)));
    tdop = sqrt(trace(Minv(4,4)));
end

function random_vecs = gen_ran_vecs(n, dim)
    rng('default');
    random_vecs = rand(n,dim) * 2 - 1;
end

function nvecs = v_normalize(vecs)
%     nvecs = vecs./vecnorm(vecs,2,2);
    for i = 1:size(vecs,1)
        nvecs(i,:) = vecs(i,:)./norm(vecs(i,:),2);
    end
end