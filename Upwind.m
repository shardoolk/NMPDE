clear all;
close all;
format long;

x0 = 0;
xN = 1;
M = 100;
N = 70;
mu = 0.001;
h = (xN - x0)/N;
delt = 2*mu*h;

for i = 1:N-1
    x(i) = i*h;
    U0(i) = sin(pi*(x(i)));
end
U1 = zeros(N-1,1);
t = 0;
for k = 1:M
    U1(1)  = U0(1) - mu*(U0(2)^2 - U0(1)^2);
    for j = 2 : N-2
        U1(j) = U0(j) - mu*(U0(j+1)^2 - U0(j)^2) ;
    end
    
    U1(N-1) =  U0(N-1) - mu*(U0(N-1)^2 - U0(N-2)^2);
    
    t = t + delt;
    U0 = U1;
end


plot(x,U1, '*')
    
    