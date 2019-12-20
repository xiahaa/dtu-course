function res = secondorder_der_gaussian_fast(t, sigma)
    x = round(-sigma*t):1:round(sigma*t);
    res = [(2^(1/2).*x.^2.*exp(-x.^2./(2*t^2)))./(2*t^4*pi^(1/2)*(t^2)^(1/2)) - ...
        (2^(1/2).*exp(-x.^2./(2*t^2)))./(2*t^2*pi^(1/2)*(t^2)^(1/2))];
end
