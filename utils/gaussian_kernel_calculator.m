function g = gaussian_kernel_calculator(D, t, sigma)
% automatically compute the Gaussian kernel by providing the dimension and
% the standard variance.
% Derivation:
% assume multivariables are isometric and independent, which means their
% covariance matrix is formed as s^2.I. Then the multivariable Gaussian
% distribution is given with the full equation as:
% 1/((2pi)^D*det(t^2.*I))^(1/2) exp(-0.5*(x-u)'*1/(t^2)*inv(I)*(x-u)) 
% simplify this by det(t^2.*I)=(t^2)^D=t^(2D) and u=0, we have
% 1/((2pi)^D*t^(2D))^(1/2) exp(-0.5*(x)'*1/(t^2)*(x)) = 
% 1/((2pi)^D*t^(2D))^(1/2) exp(-0.5*1/(t^2).*(x).^2)
% Author: xiahaa@space.dtu.dk
    if D == 1
        x = round(-t):round(t);
        f = @(x) (1/(((2*pi)^D*sigma^(2*D))^(0.5)).*exp((-0.5/(sigma^2)).*(x.^2)));
        g = f(x);
    elseif D == 2
        u = round(-t):round(t);
        [x,y] = meshgrid(u,u);
        f = @(x,y) (1/(((2*pi)^D*sigma^(2*D))^(0.5)).*exp((-0.5/(sigma^2)).*(x.^2+y.^2)));
        g = f(x,y);
    elseif D == 3
        u = round(-t):round(t);
        [x,y,z] = meshgrid(u,u,u);
        f = @(x,y,z) (1/(((2*pi)^D*sigma^(2*D))^(0.5)).*exp((-0.5/(sigma^2)).*(x.^2+y.^2+z.^2)));
        g = f(x,y,z);
    end
    g = g./sum(vec(g));
end