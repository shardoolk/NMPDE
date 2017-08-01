clear all;
close all;

format long;

a = 0;
b = 1;
t0 = 0;
tf = 1;

N = 10;
M = 1000;

del_t = (tf-t0)/M;

for p=1:5
h(p) = (b-a)/N;

mu = del_t/(2*h(p)^2);

u0 = zeros(N-1,1);
x = zeros(N-1,1);
for j=1:N-1
    x(j) = j*h(p);
    u0(j) = sin(pi*x(j));
end

F = zeros(N-1,1);
A = zeros(N-1,N-1);

U1 = zeros(N-1,1);
U = ones(N-1,1);
%%
% $x^2+e^{\pi i}$
error = 1000;
tol = 10^-4;

t=t0;
while error > tol
    for k=1:M
        F(1) = (1 + 2*mu + delt)*U0(1) - delt*(U0(1)^2) - mu*U0(2) 