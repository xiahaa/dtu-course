function varargout = navSolver(prs, sat_pos, options)
%% GNSS navigation solver: calculate antenna's position from at least 4
%% pseudoranges.
    global cspd;
    cspd = 299792458;% m / s;
    
    dim = 4;
    x0 = zeros(dim,1);
    
    %% estimate initial value
    if options.useSOCP == 1
        x0 = SOCP(prs,sat_pos);
    elseif options.useDLT == 1
        x0 = DLT(prs,sat_pos);
    elseif options.useBancroft == 1
        x0 = bancroft_fast(prs,sat_pos,[0;0;0;0]);
    elseif options.usePrior == 1
        x0 = [options.x0_prior(1:3);options.x0_prior(4)] ;%% plus 1km as requested
    else
        error('Choose either SOCP or DLT for initial estimation.');
        return;
    end
    
    if options.useWLS == 1
        [x,std_x,QDOP,Qenu,llh] = solveWLS(prs, sat_pos, x0, options);
    elseif options.useGN == 1
        [x,std_x,QDOP,Qenu,llh] = solveGaussNewton(prs, sat_pos, x0, options);
    elseif options.useSD == 1
        [x,std_x,QDOP,Qenu,llh] = solveBySteepestDescentMethod(prs, sat_pos, x0, options);
    elseif options.useLM == 1
        [x,std_x,QDOP,Qenu,llh] = solveByLevenbergMarquardt(prs, sat_pos, x0, options);
    end
    
    x(4) = x(4)./cspd;
    varargout{1} = x;
    varargout{2} = std_x;
    varargout{3} = QDOP;
    varargout{4} = Qenu;
    varargout{5} = llh;
end




function x0 = DLT(prs, sat_pos)
%% use DLT to estimate an initial value
    %% initially, omit receiver clock error
    n = size(prs,1);
    v = 1:n;
    C = nchoosek(v,2);
    cspd = 299792458;% m / s;
    sat_pos_r = sat_pos(:,1).^2+sat_pos(:,2).^2+sat_pos(:,3).^2;
    
    dprs = prs(C(:,1)).^2 - prs(C(:,2)).^2 - sat_pos_r(C(:,1),1) + sat_pos_r(C(:,2),1);
    A = 2.*[sat_pos(C(:,2),1)-sat_pos(C(:,1),1), ...
            sat_pos(C(:,2),2)-sat_pos(C(:,1),2), ...
            sat_pos(C(:,2),3)-sat_pos(C(:,1),3)];
    x0 = lssolver(A,dprs);
%     x0 = [x0;0];
    rho = sqrt((sat_pos(:,1)-x0(1)).^2+(sat_pos(:,2)-x0(2)).^2+(sat_pos(:,3)-x0(3)).^2);
    dt = mean(prs-rho) / cspd;
    x0 = [x0;dt];
    
end

function x0 = SOCP(prs,sat_pos)
%% estimate initial value using SOCP and relevant relaxation    
    sat_pos = sat_pos ./ 1000;
    prs = prs ./ 1000;

    n1 = size(prs,1);
    dim = 4;
    
    import mosek.fusion.*
    M = Model('cqo1');
    c = M.variable('c', n1, Domain.greaterThan(0.0));
    x = M.variable('x', dim, Domain.unbounded());
    t = M.variable('t', 1, Domain.greaterThan(0.0));
    
    % create the aliases
    a1 = Expr.sub(Var.repeat(x.index(1),n1),sat_pos(:,1));
    b1 = Expr.sub(Var.repeat(x.index(2),n1),sat_pos(:,2));
    c1 = Expr.sub(Var.repeat(x.index(3),n1),sat_pos(:,3));
    e = Expr.sub(prs,c);
    d = Expr.sub(e, Var.repeat(x.index(4),n1));
    %% socp cone 
    M.constraint('c1',Expr.hstack(c,Expr.hstack(a1,b1,c1)),Domain.inQCone());
    M.constraint('c2',Expr.hstack(t,Expr.constTerm(1, 1.0/2.0),Expr.transpose(d)),Domain.inRotatedQCone());
    
    regulization = 1;
    v1 = regulization.*ones(n1,1);%% in order to control the relaxation, could add regulization on c1
    M.objective("obj", ObjectiveSense.Minimize, Expr.add(t,Expr.dot(v1,c)));
    M.solve();
    
    disp(toString(M.getPrimalSolutionStatus));
    disp(num2str(t.level()));
    
    x0 = x.level();
    
    x0 = x0 .* 1000;
end

