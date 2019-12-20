
function varargout = fit_ellipse(points)
    x = points(1,:)';
    y = points(2,:)';
    n = size(x,1);
    D = [x.^2 x.*y y.^2 x y ones(n,1)];
    S = D'*D;
    C(6,6)=0;C(1,3)=-2;C(2,2)=1;C(3,1)=-2;
    [gvec,geval] = eig(S,C);
    [~,nc] = find(geval<0 & ~isinf(geval));
    params = gvec(:,nc);
    fitting_err = (D*params);
    fitting_err = fitting_err'*fitting_err;
    varargout{1} = params;
    varargout{2} = fitting_err/n;
end