function [polyCoeffs,realt]=cvx_project_traj_gen_solver_sdp_mo(p,initv,inita,endv,enda,maxv,maxa,pz,initvz,endvz,initaz,endaz,maxvz,maxaz,l)
    tic
    dims = 3;
    N = size(p,1);
    pWpts = zeros(N,3);
    vWpts = zeros(N,3);
    aWpts = zeros(N,3);
     
    pWpts(:,1:2) = p;
    vWpts(1,1) = initv(1);
    vWpts(1,2) = initv(2);
    vWpts(end,1) = endv(1);
    vWpts(end,2) = endv(2);
    vWpts(2:end-1,1:2) = inf;

    aWpts(1,1) = inita(1);
    aWpts(1,2) = inita(2);
    aWpts(end,1) = enda(1);
    aWpts(end,2) = enda(2);
    aWpts(2:end-1,1:2) = inf;
    
    pWpts(:,3) = pz;
    vWpts(1,3) = initvz;
    vWpts(end,3)=endvz;
    aWpts(1,3) = initaz;
    aWpts(end,3)=endaz;
    vWpts(2:end-1,3) = inf;
    aWpts(2:end-1,3) = inf;
    
    vzmax = maxvz;%
    azmax = maxaz;%
    
    %% compute hyperplane    
    p31 = pWpts(1:end-1,:);
    p32 = pWpts(2:end,:);
    pdir = p32-p31;
    pdir = pdir./vecnorm(pdir,2);
    n1 = [zeros(size(pdir,1),1) -pdir(:,3) pdir(:,2)];
    n1 = n1./vecnorm(n1,2,2);
    d1 = -(n1(:,1).*pWpts(1:end-1,1)+n1(:,2).*pWpts(1:end-1,2)+n1(:,3).*pWpts(1:end-1,3));
    n2 = cross(n1,pdir,2);
    n2 = n2./vecnorm(n2,2,2);
    d2 = -(n2(:,1).*pWpts(1:end-1,1)+n2(:,2).*pWpts(1:end-1,2)+n2(:,3).*pWpts(1:end-1,3));
    hyperplanes = [n1 d1 n2 d2];
    
    vmax = maxv;% 5m/s
    amax = maxa;% 6.8 m/s^2
    
    t = linspace(0,size(pWpts,1)-1,size(pWpts,1));
    order = 5;
 
    %% 1. pre-allocated timestamp, may not be reasonable
    [polyCoeffs] = trajectoryGenerator(pWpts, vWpts, aWpts, t, dims, order, l, hyperplanes,[],[]);
   
%     visulization(pWpts, polyCoeffs, t, order, vmax, amax);    
    corridorExtremes = cell(size(pWpts,1)-1,1);

    origS = [1;1];
    realt = t;
    while 1
        segments = size(pWpts,1)-1;
        numCoeff = order + 1;
        replan = 0;
        %% 2. find extreme point
        vxext = {};vyext = {};vzext = {};
        axext = {};ayext = {};azext = {};
        p1maxmin = [];p2maxmin = [];
        
        for i=1:segments
            scalar = 1./(realt(i+1)-realt(i));
            
            polyx = polyCoeffs((i-1)*numCoeff+1:i*numCoeff,1);
            polyy = polyCoeffs((i-1)*numCoeff+1:i*numCoeff,2);
            polyz = polyCoeffs((i-1)*numCoeff+1:i*numCoeff,3);

            %% v
            vxpoly = scalar.*[5*polyx(6) 4*polyx(5) 3*polyx(4) 2*polyx(3) polyx(2)];
            vypoly = scalar.*[5*polyy(6) 4*polyy(5) 3*polyy(4) 2*polyy(3) polyy(2)];
            vzpoly = scalar.*[5*polyz(6) 4*polyz(5) 3*polyz(4) 2*polyz(3) polyz(2)];
            %% a
            axpoly = scalar^2.*[20*polyx(6) 12*polyx(5) 6*polyx(4) 2*polyx(3)];
            aypoly = scalar^2.*[20*polyy(6) 12*polyy(5) 6*polyy(4) 2*polyy(3)];
            azpoly = scalar^2.*[20*polyz(6) 12*polyz(5) 6*polyz(4) 2*polyz(3)];
            %% vdot
            vdxpoly = axpoly;
            vdypoly = aypoly;
            vdzpoly = azpoly;
            %% adot
            adxpoly = scalar^3.*[60*polyx(6) 24*polyx(5) 6*polyx(4)];
            adypoly = scalar^3.*[60*polyy(6) 24*polyy(5) 6*polyy(4)];
            adzpoly = scalar^3.*[60*polyz(6) 24*polyz(5) 6*polyz(4)];
            %% roots of vdot
            xroots = roots(vdxpoly);
            yroots = roots(vdypoly);
            zroots = roots(vdzpoly);
            %% validation
            idx = isreal(xroots) & real(xroots) > 0 & real(xroots) < 1;
            xroots = xroots(idx);
            idy = isreal(yroots) & real(yroots) > 0 & real(yroots) < 1;
            yroots = yroots(idy);
            idz = isreal(zroots) & real(zroots) > 0 & real(zroots) < 1;
            zroots = zroots(idz);
            %% possible extreme points
            xroots = [xroots;0;1];
            yroots = [yroots;0;1];
            zroots = [zroots;0;1];
            %% current value
            if ~isempty(xroots)
                vxmaxmin = polyval(vxpoly,xroots);
            end
            if ~isempty(yroots)
                vymaxmin = polyval(vypoly,yroots);
            end
            if ~isempty(zroots)
                vzmaxmin = polyval(vzpoly,zroots);
            end
            %% roots of ador
            xroots = roots(adxpoly);
            yroots = roots(adypoly);
            zroots = roots(adzpoly);
            %% validataion
            idx = isreal(xroots) & real(xroots) > 0 & real(xroots) < 1;
            xroots = xroots(idx);
            idy = isreal(yroots) & real(yroots) > 0 & real(yroots) < 1;
            yroots = yroots(idy);
            idz = isreal(zroots) & real(zroots) > 0 & real(zroots) < 1;
            zroots = zroots(idz);
            %% possible extreme points
            xroots = [xroots;0;1];
            yroots = [yroots;0;1];
            zroots = [zroots;0;1];
            %% current value
            if ~isempty(xroots)
                axmaxmin = polyval(axpoly,xroots);
            end
            if ~isempty(yroots)
                aymaxmin = polyval(aypoly,yroots);
            end
            if ~isempty(zroots)
                azmaxmin = polyval(azpoly,zroots);
            end
            %% store
            vxext{i} = vxmaxmin;
            vyext{i} = vymaxmin;
            vzext{i} = vzmaxmin;
            %% store
            axext{i} = axmaxmin;
            ayext{i} = aymaxmin;
            azext{i} = azmaxmin;
            %% if already satisfactory
            r1 = isempty(find(abs(vxmaxmin) > vmax)) && ...
                     isempty(find(abs(vymaxmin) > vmax)) && ...
                     isempty(find(abs(axmaxmin) > amax)) && ...
                     isempty(find(abs(aymaxmin) > amax)) && ...
                     isempty(find(abs(vzmaxmin) > vzmax)) && ...
                     isempty(find(abs(azmaxmin) > azmax));
            replan = replan + r1;
        end

        lambda = 10;
        
        if replan == segments
            break;
        end
        
        %% try to find an linear mapping scalar for each segment, s1, s2, ..., sn
        %% the scalar is used to map from real world timestamp to virtual timestamp
        scalars = [];
        
        %% USE CVX, model this as a SDP problem with SDR
        %% firstly reformulate cost function
        Aopt = [];
        bopt = [];
        nopt = 4;

        %% vx
        for j = 1:segments
            polyx = polyCoeffs((j-1)*numCoeff+1:j*numCoeff,1);
            polyy = polyCoeffs((j-1)*numCoeff+1:j*numCoeff,2);
            polyz = polyCoeffs((j-1)*numCoeff+1:j*numCoeff,3);

            C = costC(polyx, polyy, polyz);
            C0 = C+[1 -2 1 0;0 0 0 0;0 0 0 0;0 0 0 0].*lambda;            

            kk = 1;

            [Aopt,bopt,kk] = addAugumentedConstraint(vxext,Aopt,bopt,kk,j,vmax,'v');
            [Aopt,bopt,kk] = addAugumentedConstraint(vyext,Aopt,bopt,kk,j,vmax,'v');
            [Aopt,bopt,kk] = addAugumentedConstraint(vzext,Aopt,bopt,kk,j,vzmax,'v');

            [Aopt,bopt,kk] = addAugumentedConstraint(axext,Aopt,bopt,kk,j,amax,'a');
            [Aopt,bopt,kk] = addAugumentedConstraint(ayext,Aopt,bopt,kk,j,amax,'a');
            [Aopt,bopt,kk] = addAugumentedConstraint(azext,Aopt,bopt,kk,j,azmax,'a');

            cvx_begin quiet
                variable X(nopt, nopt) hankel
                minimize( trace(C0*X) );
                subject to 
                    for ii = 1:kk-1
                        trace(Aopt(:,:,ii)*X) <= bopt(ii);
                    end
                    X >= 0
                    X(1,1) == 1
                    X(1,2) >= 0
                    X(1,3) >= 0
                    X(1,4) >= 0
                    X(2,4) >= 0
                    X(3,4) >= 0
                    X(4,4) >= 0
            cvx_end

            %% here, I just take the maximum eigvalue corresponding solution
            [eigvec, eigval] = eig(X);
            maxEig = eigval(end,end);
            correspondingEvec = eigvec(:,end);
            xsol = sqrt(maxEig).*correspondingEvec;
            %% parsing solution
            tt = double(xsol(2)/xsol(1));
            scalars = [scalars;tt];
        end
        
        ttmp = realt;
        for i = 2:numel(t)
            realt(i) = realt(i-1) + (ttmp(i)-ttmp(i-1)) / scalars(i-1);
        end
        origS = scalars;
        [polyCoeffs] = trajectoryGenerator(pWpts, vWpts, aWpts, realt, dims, order, l, hyperplanes,[],[]);
    end
%     visulization(pWpts, polyCoeffs, realt, order, vmax, amax);
    toc
end

function [Aopt,bopt,kk] = addAugumentedConstraint(vext,Aopt,bopt,kk,j,vmax,derv)

    if derv == 'v'
        for k = 1:size(vext{j},1)
            Aopt(:,:,kk) = [0 vext{j}(k) 0 0; 0 0 0 0;0 0 0 0;0 0 0 0];
            bopt(kk) = vmax;
            Aopt(:,:,kk+1) = [0 -vext{j}(k) 0 0; 0 0 0 0;0 0 0 0;0 0 0 0];
            bopt(kk+1) = vmax;
            kk = kk + 2;
                    
            % augumented constraints
            Aopt(:,:,kk) = [0 -vmax vext{j}(k) 0; 0 0 0 0;0 0 0 0;0 0 0 0];
            bopt(kk) = 0;
            Aopt(:,:,kk+1) = [0 -vmax -vext{j}(k) 0; 0 0 0 0;0 0 0 0;0 0 0 0];
            bopt(kk+1) = 0;
            kk = kk + 2;
                    
            Aopt(:,:,kk) = [0 0 -vmax vext{j}(k); 0 0 0 0;0 0 0 0;0 0 0 0];
            bopt(kk) = 0;
            Aopt(:,:,kk+1) = [0 0 -vmax -vext{j}(k); 0 0 0 0;0 0 0 0;0 0 0 0];
            bopt(kk+1) = 0;
            kk = kk + 2;
                    
            Aopt(:,:,kk) = [0 0 0 -vmax; 0 0 0  vext{j}(k);0 0 0 0;0 0 0 0];
            bopt(kk) = 0;
            Aopt(:,:,kk+1) = [0 0 0 -vmax; 0 0 0  -vext{j}(k);0 0 0 0;0 0 0 0];
            bopt(kk+1) = 0;
            kk = kk + 2;
                    
            Aopt(:,:,kk) = [0 0 0 0; 0 0 0 -vmax;0 0 0 vext{j}(k);0 0 0 0];
            bopt(kk) = 0;
            Aopt(:,:,kk+1) = [0 0 0 0; 0 0 0 -vmax;0 0 0 -vext{j}(k);0 0 0 0];
            bopt(kk+1) = 0;
            kk = kk + 2;
                    
            Aopt(:,:,kk) = [0 0 0 0; 0 0 0 0;0 0 0 -vmax;0 0 0 vext{j}(k)];
            bopt(kk) = 0;
            Aopt(:,:,kk+1) = [0 0 0 0; 0 0 0 0;0 0 0 -vmax;0 0 0 -vext{j}(k)];
            bopt(kk+1) = 0;
            kk = kk + 2;
        end
    elseif derv == 'a'
        for k = 1:size(vext{j},1)
            Aopt(:,:,kk) = [0 0 vext{j}(k) 0;0 0 0 0;0 0 0 0;0 0 0 0];
            bopt(kk) = vmax;
            Aopt(:,:,kk+1) = [0 0 -vext{j}(k) 0;0 0 0 0;0 0 0 0;0 0 0 0];
            bopt(kk+1) = vmax;
            kk = kk + 2;

            Aopt(:,:,kk) = [0 -vmax 0 vext{j}(k);0 0 0 0;0 0 0 0;0 0 0 0];
            bopt(kk) = 0;
            Aopt(:,:,kk+1) = [0 -vmax 0 -vext{j}(k);0 0 0 0;0 0 0 0;0 0 0 0];
            bopt(kk+1) = 0;
            kk = kk + 2;

            Aopt(:,:,kk) = [0 0 -vmax 0;0 0 0 vext{j}(k);0 0 0 0;0 0 0 0];
            bopt(kk) = 0;
            Aopt(:,:,kk+1) = [0 0 -vmax 0;0 0 0 -vext{j}(k);0 0 0 0;0 0 0 0];
            bopt(kk+1) = 0;
            kk = kk + 2;

            Aopt(:,:,kk) = [0 0 0 -vmax;0 0 0 0;0 0 0 vext{j}(k);0 0 0 0];
            bopt(kk) = 0;
            Aopt(:,:,kk+1) = [0 0 0 -vmax;0 0 0 0;0 0 0 -vext{j}(k);0 0 0 0];
            bopt(kk+1) = 0;
            kk = kk + 2;

            Aopt(:,:,kk) = [0 0 0 0;0 0 0 -vmax;0 0 0 0;0 0 0 vext{j}(k)];
            bopt(kk) = 0;
            Aopt(:,:,kk+1) = [0 0 0 0;0 0 0 -vmax;0 0 0 0;0 0 0 -vext{j}(k)];
            bopt(kk+1) = 0;
            kk = kk + 2;
        end
    end
end

function [polyCoeffs] = trajectoryGenerator(pWpts, vWpts, aWpts, t, dims, order, l, hyperplanes, vmax, amax)
    %% order of trajectory 5
    numCoeff = order + 1;
    segments = size(pWpts,1)-1;
    polyCoeffs = zeros(segments*numCoeff,dims);

    %% cost matrix
    Q = zeros(numCoeff*segments, numCoeff*segments);
    
    %% equality constraint
    maxNumCons = numCoeff * segments;
    alreadySetCons = 0;
    
    cc1 = [1;2.*ones(size(pWpts,1)-2,1);1];
    cc2 = [1.*ones(size(pWpts,1)-2,1)];
    
    pc1 = sum(double(pWpts(:,1)~=inf).*cc1);
    vc1 = sum(double(vWpts(:,1)~=inf).*cc1);
    ac1 = sum(double(aWpts(:,1)~=inf).*cc1);
    
    ctc1 = sum(double(vWpts(2:end-1,1)==inf).*cc2);
    ctc2 = sum(double(aWpts(2:end-1,1)==inf).*cc2);
    alreadySetCons = pc1 + vc1 + ac1 + ctc1 + ctc2;
    
    %% cnstruct A
    Aeq = zeros(alreadySetCons,numCoeff*segments);
    beqs = zeros(alreadySetCons,dims);
    k = 1;
    for i=1:segments
        t0 = 0;
        t1 = 1;%t(i+1)-t(i);
        scalar = 1./(t(i+1)-t(i));
        %% start point
        if pWpts(i,1) ~= inf
            Aeq(k,(i-1)*numCoeff+1:i*numCoeff) = [1 0 0 0 0 0];
            beqs(k,:) = pWpts(i,:);
            k = k + 1;
        end
        if vWpts(i,1) ~= inf
            Aeq(k,(i-1)*numCoeff+1:i*numCoeff) = scalar.*[0 1 0 0 0 0];
            beqs(k,:) = vWpts(i,:);
            k = k + 1;
        end
        if aWpts(i,1) ~= inf
            Aeq(k,(i-1)*numCoeff+1:i*numCoeff) = scalar^2.*[0 0 2 0 0 0];
            beqs(k,:) = aWpts(i,:);
            k = k + 1;
        end
        %% end point
        if pWpts(i+1,1) ~= inf
            Aeq(k,(i-1)*numCoeff+1:i*numCoeff) = [1 1 1 1 1 1];
            beqs(k,:) = pWpts(i+1,:);
            k = k + 1;
        end
        if vWpts(i+1,1) ~= inf
            Aeq(k,(i-1)*numCoeff+1:i*numCoeff) = scalar.*[0 1 2 3 4 5];
            beqs(k,:) = vWpts(i+1,:);
            k = k + 1;
        end
        if aWpts(i+1,1) ~= inf
            Aeq(k,(i-1)*numCoeff+1:i*numCoeff) = scalar^2.*[0 0 2 6 12 20];
            beqs(k,:) = aWpts(i+1,:);
            k = k + 1;
        end
    end
    %% continuity
    for i=2:size(pWpts,1)-1
        t0 = 1;
        t1 = 0;
        s1 = 1./(t(i)-t(i-1));
        s2 = 1./(t(i+1)-t(i));
        if vWpts(i,1) == inf
            Aeq(k,(i-2)*numCoeff+1:(i-1)*numCoeff) = s1.*[0 1 2 3 4 5];
            Aeq(k,(i-1)*numCoeff+1:(i)*numCoeff) = s2.*[-0 -1 0 0 0 0];
            beqs(k,:) = zeros(1,dims);
            k = k + 1;
        end
        if aWpts(i,1) == inf
            Aeq(k,(i-2)*numCoeff+1:(i-1)*numCoeff) = s1^2.*[0 0 2 6 12 20];
            Aeq(k,(i-1)*numCoeff+1:(i)*numCoeff) = s2^2.*[-0 -0 -2 0 0 0];
            beqs(k,:) = zeros(1,dims);
            k = k + 1;
        end
    end
    
    %% Q
    for i = 1:segments
        Qs = constructQ(t(i+1)-t(i));
        Q((i-1)*numCoeff+1:(i)*numCoeff, (i-1)*numCoeff+1:(i)*numCoeff) = Qs;
    end

    %% prepare inequality constraints
    Aineq = [];
    bineq = [];
    nn = 5;
    for i = 1:segments
        hyperplane = hyperplanes(i,:);
        ts = linspace(0, 1, nn)';
        numIneq = size(numel(ts),2);
        Ax = zeros(nn,numCoeff*segments);
        Ay = zeros(nn,numCoeff*segments);
        Az = zeros(nn,numCoeff*segments);
        
        %% position constraint
        ttp = [ones(nn,1) ts ts.^2 ts.^3 ts.^4 ts.^5];
        Ax(:,(i-1)*numCoeff+1:i*numCoeff) = hyperplane(1).*ttp;
        Ay(:,(i-1)*numCoeff+1:i*numCoeff) = hyperplane(2).*ttp;
        Az(:,(i-1)*numCoeff+1:i*numCoeff) = hyperplane(3).*ttp;
        
        Aineq = [Aineq;[Ax Ay Az];[-Ax -Ay -Az]];
        bineq = [bineq;repmat(0.5*l-hyperplane(4),nn,1);repmat(0.5*l+hyperplane(4),nn,1)];
        
        Ax(:,(i-1)*numCoeff+1:i*numCoeff) = hyperplane(5).*ttp;
        Ay(:,(i-1)*numCoeff+1:i*numCoeff) = hyperplane(6).*ttp;
        Az(:,(i-1)*numCoeff+1:i*numCoeff) = hyperplane(7).*ttp;
        
        Aineq = [Aineq;[Ax Ay Az];[-Ax -Ay -Az]];
        bineq = [bineq;repmat(0.5*l-hyperplane(8),nn,1);repmat(0.5*l+hyperplane(8),nn,1)];

        if ~isempty(vmax)
            %% velocity constaint
            Ax = zeros(nn,numCoeff*segments);
            Ay = zeros(nn,numCoeff*segments);
            Az = zeros(nn,numCoeff*segments);

            s1 = 1./(t(i+1)-t(i));
            ttv = s1.*[zeros(nn,1) ones(nn,1) 2.*ts 3.*ts.^2 4.*ts.^3 5.*ts.^4];
            Ax(:,(i-1)*numCoeff+1:i*numCoeff) = ttv;
            Ay(:,(i-1)*numCoeff+1:i*numCoeff) = ttv;
            Az(:,(i-1)*numCoeff+1:i*numCoeff) = ttv;
            Aineq = [Aineq;[Ax Ay Az];[-Ax -Ay -Az]];
            bineq = [bineq;repmat(vmax,2*nn,1)];
        end
        if ~isempty(amax)
            Ax = zeros(nn,numCoeff*segments);
            Ay = zeros(nn,numCoeff*segments);
            Az = zeros(nn,numCoeff*segments);

            tta = (s1*s1).*[zeros(nn,1) zeros(nn,1) 2.*ones(nn,1) 6.*ts 12.*ts.^2 20.*ts.^3];
            Ax(:,(i-1)*numCoeff+1:i*numCoeff) = tta;
            Ay(:,(i-1)*numCoeff+1:i*numCoeff) = tta;
            Az(:,(i-1)*numCoeff+1:i*numCoeff) = tta;
            Aineq = [Aineq;[Ax Ay Az];[-Ax -Ay -Az]];
            bineq = [bineq;repmat(amax,2*nn,1)];
        end
    end
    Q1 = blkdiag(Q,Q,Q);
    Aeq1 = blkdiag(Aeq,Aeq,Aeq);
    beq1 = [beqs(:,1);beqs(:,2);beqs(:,3)];

    sols = quadprog(2*Q1,[],Aineq,bineq,Aeq1,beq1);

    xsol = sols(1:size(Q,1));
    ysol = sols(size(Q,1)+1:size(Q,1)*2);
    zsol = sols(size(Q,1)*2+1:size(Q,1)*3);

    polyCoeffs(:,1) = xsol(1:size(Q,1));
    polyCoeffs(:,2) = ysol(1:size(Q,1));
    polyCoeffs(:,3) = zsol(1:size(Q,1));
end

function p = costPoly(poly, x, ts)
    p3 = poly(4);
    p4 = poly(5);
    p5 = poly(6);

    a1 = 36*p3^2;
    a2 = 144*p3*p4;
    a3 = (192*p4^2 + 240*p3*p5);
    a4 = 720*p4*p5;
    a5 = 720*p5^2;
    
    p = a1*x+a2*x^2+a3*x^3+a4*x^4+a5*x^5 - (a1*ts+a2*ts^2+a3*ts^3+a4*ts^4+a5*ts^5);
end

function sol = solveViaQR(Q,A,b)
    [Qt Rt] = qr(A');
    p = size(A,1);
    n = size(A,2);
    Q1 = Qt(1:n,1:p);
    Q2 = Qt(1:n,p+1:end);
    R = Rt(1:p,1:p);
    x0 = A\((R'*Q1'*Q1*inv(R)')*b);
    Q11 = Q2'*Q*Q2;
    Q12 = Q*Q2;
    z = -(Q11)\(Q12'*x0);
    sol = x0 + Q2*z;
end

function Qs = constructQ(dt)
    scalar = 1/dt;
    Qs = scalar^3.*[ 0, 0, 0,   0,   0,   0; ...
                       0, 0, 0,   0,   0,   0; ...
                       0, 0, 0,   0,   0,   0; ...
                       0, 0, 0,  36,  72, 120; ...
                       0, 0, 0,  72, 192, 360; ...
                       0, 0, 0, 120, 360, 720];
end

function C = costC(polyx, polyy, polyz)
    p3 = polyx(4);
    p4 = polyx(5);
    p5 = polyx(6);

    a1 = 36*p3^2;
    a2 = 144*p3*p4;
    a3 = (192*p4^2 + 240*p3*p5);
    a4 = 720*p4*p5;
    a5 = 720*p5^2;
    
    ps = a1+a2+a3+a4+a5;
    
    Cx = [0 0 0 0;...
         0 0 0 0;...
         0 0 0 0;...
         0 0 ps 0];
     
    p3 = polyy(4);
    p4 = polyy(5);
    p5 = polyy(6);

    a1 = 36*p3^2;
    a2 = 144*p3*p4;
    a3 = (192*p4^2 + 240*p3*p5);
    a4 = 720*p4*p5;
    a5 = 720*p5^2;
    ps = a1+a2+a3+a4+a5;
    Cy = [0 0 0 0;...
         0 0 0 0;...
         0 0 0 0;...
         0 0 ps 0];
     
    p3 = polyz(4);
    p4 = polyz(5);
    p5 = polyz(6);

    a1 = 36*p3^2;
    a2 = 144*p3*p4;
    a3 = (192*p4^2 + 240*p3*p5);
    a4 = 720*p4*p5;
    a5 = 720*p5^2;
    ps = a1+a2+a3+a4+a5;
    Cz = [0 0 0 0;...
         0 0 0 0;...
         0 0 0 0;...
         0 0 ps 0];
     
     
    C = Cx + Cy + Cz;
end

function [tss,pts,vts,ats] = sample_pvas(pWpts, polyCoeffs, t, order, cnt)
    numCoeff = order+1;
    tss = [];
    pts = [];
    vts = [];
    ats = [];
    for i = 1:size(t,2)-1
        s = 1./(t(i+1)-t(i));
        ts = linspace(0,1,cnt)';
        vtp = [ones(numel(ts),1) ts ts.^2 ts.^3 ts.^4 ts.^5];
        vtv = [zeros(cnt,1) ones(cnt,1) 2.*ts 3.*ts.^2 4.*ts.^3 5.*ts.^4];
        vta = [zeros(cnt,1) zeros(cnt,1) 2.*ones(cnt,1) 6.*ts 12.*ts.^2 20.*ts.^3];
        
        polyx = polyCoeffs((i-1)*numCoeff+1:i*numCoeff,1);
        polyy = polyCoeffs((i-1)*numCoeff+1:i*numCoeff,2);
        polyz = polyCoeffs((i-1)*numCoeff+1:i*numCoeff,3);
        
        xsamples = vtp*polyx;
        ysamples = vtp*polyy;
        zsamples = vtp*polyz;
        
        vxsamples = vtv*polyx.*s;
        vysamples = vtv*polyy.*s;
        vzsamples = vtv*polyz.*s;
        
        axsamples = vta*polyx.*(s^2);
        aysamples = vta*polyy.*(s^2);
        azsamples = vta*polyz.*(s^2);
        
        pts = [pts;[xsamples ysamples zsamples]];
        vts = [vts;[vxsamples vysamples vzsamples]];
        ats = [ats;[axsamples aysamples azsamples]];
        
        tss = [tss;ts+(i-1)];
    end
end

function visulization(pWpts, polyCoeffs, t, order, vmax, amax)
    figure
    plot(pWpts(:,1),pWpts(:,2),'ro','MarkerSize',5);
    hold on;
    grid on;
    pts = [];
    vts = [];
    ats = [];

    numCoeff = order+1;
    tss = [];
    for i = 1:size(t,2)-1
        scalar = 1./(t(i+1)-t(i));
        ts = linspace(0,1,1000)';
        polyx = polyCoeffs((i-1)*numCoeff+1:i*numCoeff,1);
        polyy = polyCoeffs((i-1)*numCoeff+1:i*numCoeff,2);
        xsamples = polyx(1) + polyx(2).*ts + polyx(3).*ts.^2 + polyx(4).*ts.^3 + polyx(5).*ts.^4 + polyx(6).*ts.^5;
        ysamples = polyy(1) + polyy(2).*ts + polyy(3).*ts.^2 + polyy(4).*ts.^3 + polyy(5).*ts.^4 + polyy(6).*ts.^5;
        
        vxsamples = polyx(2) + polyx(3).*ts.*2 + polyx(4).*3.*ts.^2 + polyx(5).*4.*ts.^3 + polyx(6).*5.*ts.^4;
        vysamples = polyy(2) + polyy(3).*ts.*2 + polyy(4).*3.*ts.^2 + polyy(5).*4.*ts.^3 + polyy(6).*5.*ts.^4;
        
        axsamples = polyx(3).*2 + polyx(4).*6.*ts.^1 + polyx(5).*12.*ts.^2 + polyx(6).*20.*ts.^3;
        aysamples = polyy(3).*2 + polyy(4).*6.*ts.^1 + polyy(5).*12.*ts.^2 + polyy(6).*20.*ts.^3;
        
        vxsamples = vxsamples.*scalar;
        vysamples = vysamples.*scalar;
        axsamples = axsamples.*(scalar^2);
        aysamples = aysamples.*(scalar^2);
        
        pts = [pts;[xsamples ysamples]];
        vts = [vts;[vxsamples vysamples]];
        ats = [ats;[axsamples aysamples]];
        
        tss = [tss;ts+(i-1)];
    end
    
    plot(pts(:,1),pts(:,2),'b-');
    title('Trajectory');
    figure
    subplot(2,1,1);
    plot(tss,vts(:,1),'b-');grid on;hold on;
    plot(tss,ones(numel(vts(:,1)),1).*vmax,'r.-');
    plot(tss,ones(numel(vts(:,1)),1).*-vmax,'r.-');
    legend('Generated','Maximum');
    title('Velocity');
    subplot(2,1,2);
    plot(tss,vts(:,2),'b-');grid on;hold on;
    plot(tss,ones(numel(vts(:,2)),1).*vmax,'r.-');
    plot(tss,ones(numel(vts(:,2)),1).*-vmax,'r.-');
    legend('Generated','Maximum');
    figure    
    subplot(2,1,1);
    plot(tss,ats(:,1),'b-');grid on;hold on;
    plot(tss,ones(numel(ats(:,1)),1).*amax,'r.-');
    plot(tss,ones(numel(ats(:,1)),1).*-amax,'r.-');
    legend('Generated','Maximum');
    title('Acceleration');
    subplot(2,1,2);
    plot(tss,ats(:,2),'b-');grid on;hold on;
    plot(tss,ones(numel(ats(:,2)),1).*amax,'r.-');
    plot(tss,ones(numel(ats(:,2)),1).*-amax,'r.-');
    legend('Generated','Maximum');
end



