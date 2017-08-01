clear all
close all
clc
N = 10;

format long;
for p = 1:5
    x0 = 0;
    xN = 1;
     t0 = 0;
    tf = 1;
     delt = 0.0001;
M = fix((tf-t0)/delt);
   h(p) = (xN - x0)/N;

lambda = delt/(h(p)^2);
    
    for i = 1:N-1
    x(i) = i*h(p);
    U0(i) = cos(pi*(x(i) - 0.5));
end
A = zeros(N-1,N-1);
t = t0;
b = zeros(N-1,1);
for j = 1:N-1
    A(j,j) = 1 + 2*lambda;
end
for i = 1:N-1
    b(i) = U0(i);
end
for j = 1:N-2
    A(j+1,j) = -1*lambda;
    A(j,j+1) = -1*lambda;
end

for i = 1:M
    t = t + delt;
    for k = 1:N-1
        b(k) = U0(k) + delt*fcts(t,x(k));
    end 
    
    Asol = A\b;
    U0 = Asol;
end


for i = 1:N-1
    exact(i) = exp(-t)*cos(pi*(x(i) - 0.5));
end 

 




[Asol , exact']

error(p) = max(abs(exact' - Asol));

N = N*2;

end

for i = 1:p-1
    order(i) = log(error(i)/error(i+1))/log(2);
end
order

plot(x,Asol, '*', x,exact)

plot(log(h),log(error));

