function sym_compute_gradient

syms phi lambda h x y z real

f = 1/298.257223563;
e = sqrt(2*f-f^2);
a = 6378137.0;

N = a^2 / (1-e^2*sin(phi)^2)^(1/2);
p = sqrt(x^2+y^2);

phi = atan2(z, p*(1-e^2*(N/(N+h))));

f = [tan(phi);tan(lambda);h] - [z/(p*(1-e^2*(N/(N+h))));y/x;p/cos(phi)-N];
F = f'*f;

diff(F,phi)
diff(F,lambda)
diff(F,h)


end