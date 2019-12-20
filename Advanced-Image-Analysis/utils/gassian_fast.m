function res = gassian_fast(t, sigma)
    x = round(-sigma*t):1:round(sigma*t);
    res = (2^(1/2).*exp(-x.^2./(2*t^2)))./(2*pi^(1/2)*(t^2)^(1/2));
end

