function induction_time_int

syms t real
assume(t,'positive');

syms p0 p1 p2 p3 p4 p5 real

P = [p0;p1;p2;p3;p4;p5];
f0 = [1 t t^2 t^3 t^4 t^5] * P;
f1 = diff(f0,t);
f2 = diff(f1,t);
f3 = diff(f2,t);

J = int(f3*f3,t);

end