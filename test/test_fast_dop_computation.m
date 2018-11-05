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

    
    dop1 = zeros(runcnt,1);
    dop2 = zeros(runcnt,1);
    dop3 = zeros(runcnt,1);
    dop4 = zeros(runcnt,1);
    
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
        dop4(k) = eig_decomp2(H);
        t4s(k) = toc;
    end
    
    
    
    mean(t1s)
    mean(t2s)
    mean(t3s)
    mean(t4s)
    
    sum(t1s)
    sum(t2s)
    sum(t3s)
    sum(t4s)
    
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

function dop = eig_decomp(H)
    E = eig(H'*H);
    dop = sqrt(1/E(1) + 1/E(2) + 1/E(3) + 1/E(4));
end

function dop = dop_via_cl(H)
    M = H'*H;
    h1 = trace(M);
    h2 = trace(M*M);
    h3 = trace(M*M*M);
    h4 = det(M);;
    dop = sqrt((0.5*h1*h1*h1-1.5*h1*h2+h3)/(3*h4));
end

function dop = dop_via_inverse(H)
    dop = sqrt(trace(inv(H'*H)));
end

function random_vecs = gen_ran_vecs(n, dim)
    rng('shuffle');
    random_vecs = rand(n,dim) * 2 - 1;
end

function nvecs = v_normalize(vecs)
    nvecs = vecs./vecnorm(vecs,2,2);
end