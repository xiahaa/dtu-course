function plotellipse(u,v,a,b,alpha)
    % plot ellipse
    t = linspace(0,2*pi,360);
    A = [cos(alpha) sin(alpha);-sin(alpha) cos(alpha)];
    tc = [u;v];
    pts = [a*cos(t);b*sin(t)];
    ptstrans = A * pts + repmat(tc,1,size(pts,2));
%     figure(2);
    plot(ptstrans(1,:),ptstrans(2,:),'m-');
end