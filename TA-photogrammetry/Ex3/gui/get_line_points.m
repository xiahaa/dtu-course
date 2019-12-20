function pts=get_line_points(l,sz)
    a=l(1); b=l(2); c=l(3);
    h=sz(1); w=sz(2);

    % This might cause 'divide by zero' warning:
    b = b + 1e-12;
    ys=c/-b ;
    yf=-(a*w+c)/b;
    xs=c/-a;
    xf=-(b*h+c)/a;

    m1 = [[xs;1] [xf;h] [1;ys] [w;yf]];
    w2 = [(xs<=w & xs>=1) (xf<=w & xf>=1) (ys<=h & ys>=1) (yf<=h & yf>=1)];
    v = w2>0;
    pts = [m1(:,v)];
end