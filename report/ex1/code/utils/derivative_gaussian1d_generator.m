function dg = derivative_gaussian1d_generator(t, sigma, order)
    cc = 1/sqrt(2*pi*t*t);
    syms x real;
    fg = cc.*exp((-0.5/(t*t)*x^2));
    dfg = diff(fg,order);
    xreal = round(-sigma*t):1:round(sigma*t);
    dg = double(subs(dfg, xreal));
end