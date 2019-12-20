function test_trajection_generation
clear all;
clc;

syms t real 
assume(t,'positive')

% Q = [0 0 0 0 0 0; ...
%      0 0 0 0 0 0; ...
%      0 0 4*t 6*t^2 8*t^3 10*t^4; ...
%      0 0 6*t^2 12*t^3 72/4*t^4 120/5*t^5; ...
%      0 0 8*t^3 72/4*t^4 144/5*t^5 240/6*t^6; ...
%      0 0 10*t^4 120/5*t^5 240/6*t^6 400/7*t^7];
%  
% eigv = eig(Q);

f3d = [0 0 0 6 24*t 60*t^2];
f4d = [0 0 0 0 24 120*t];

Q1 = f3d'*f3d;
Q = int(Q1,t);
% Qs = subs(Q2,ts);

A = [1 0 0 0 0 0;...
     0 1 0 0 0 0;...
     0 0 2 0 0 0;...
     1 t t^2 t^3 t^4 t^5];
%      0 1 2*t 3*t*t 4*t*t*t 5*t*t*t*t];
%      0 0 2 6*t 12*t*t 24*t*t*t];
 
QT = [Q A';A zeros(size(Q,1)+size(A,1)-size(A,2))];

end