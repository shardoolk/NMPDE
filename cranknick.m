clear all
close all
clc
format long

x0 = 0;
xN = 1;
N=10;
    
for p = 1:5
t0 = 0;
tf = 1;
delt = 0.00001;
h(p) = (xN-x0)/N;

nt = fix((tf-t0)/delt);
U0 = zeros(N+1,1);
U1 = zeros(N+1,1);
k = zeros(N+1,1);
A = zeros(N+1,N+1);
B = zeros(N+1,N+1);

tau = delt/(h(p)^2);
for i = 1:N+1
    x(i) = (i-1)*h(p);
    U0(i) = cos(pi*x(i));
end

A(1,1) = 1 + tau;
A(1,2) = -tau;
A(2,1) = -tau/2;
for j = 2:N
    A(j,j) = 1+ tau;
    A(j,j+1)= -tau/2;
    A(j+1,j) = -tau/2;
end
A(N+1,N) = -tau;
A(N+1,N+1) = 1 + tau;


B(1,1) = 1 - tau;
B(1,2) = tau;
B(2,1) = tau/2;


for j = 2:N
    B(j,j) = 1- tau;
    B(j,j+1)= tau/2;
    B(j+1,j) = tau/2;
end
B(N+1,N) = tau;
B(N+1,N+1) = 1 - tau;

t = 0;
for j = 0:nt
    for i = 1:N+1
        k(i) = delt*exp(-(delt*(j+0.5)))*cos(pi*x(i))*(pi^2-1);
    end
    Rh = B*U0 + k;
    U1 = A\Rh;
    U0 = U1;
    t = t + delt;
end
 exact = exp(-1)*cos(pi*x);
error(p) = max(abs(exact - U1'));
N = 2*N;
end
for p=2:5
   order(p) = log(error(p-1)/error(p))/log(2);
end

order(2:5)
plot(x,U1,'*' ,x, exact)
![U1  exact']
max(error)