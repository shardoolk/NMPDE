%code for ADI method

close all
clear all
clc
a = 0;
b = 1;
c = 0;
d = 1;
N = 10;
M = N;
h = (b-a)/N;
tf = 1;
x = a:h:b;
y = c:h:d; 


dt = 0.001;
mu = dt/(2*h^2);

for i = 1:N+1
    for j = 1:M+1
        u1(i,j) = sin(pi*x(i))*sin(pi*y(j));
    end
end
Tf = fix(tf/dt);

for k = 1:Tf

    for i = 1:M+1
        u2(i,1) = 0;
        u2(i,N+1) = 0;
    end
    A = sparse(M-1,M-1); f1 = zeros(M-1,1);
    
 for j = 2:N  %y sweep
        for i = 2:M
            f1(i-1) = mu*(u1(i,j-1) - 2*u1(i,j) + u1(i,j+1)) + u1(i,j);
        end
 for i = 1:M-2
     A(i,i) = 1 + 2*mu;
     A(i,i+1) = -mu;
     A(i+1,i) = -mu;
 end
 A(M-1,M-1) = 1 + 2*mu;
 
 ut = A\f1;
 for i = 1:M-1
     u2(i+1,j) = ut(i);
 end
 
 end
 
 for j = 1:N+1
     u1(1,i) = 0;
     u1(M+1,i) = 0;
 end
 
 B = sparse(N-1,N-1); f2 = zeros(N-1,1);
 
 for i = 2:N
     for j = 2:M
         f2(j-1) = mu*(u2(i-1,j) - 2*u2(i,j) + u2(i+1,j)) + u2(i,j);
     end
  for i = 1:N-2
     B(i,i) = 1 + 2*mu;
     B(i,i+1) = -mu;
     B(i+1,i) = -mu;
 end
 B(N-1,N-1) = 1 + 2*mu;
 
 ut = B\f2;
 
 
for j = 1:N-1
    u1(i,j+1) = ut(j);
end
 end
% t = t+dt;
end
for i=1:M+1,
    for j=1:N+1,
       ue(i,j) = exp(-2*pi^2*tf)*sin(pi*x(i))*sin(pi*y(j)); 
  end
  end
  
  error = max(max(abs(u1-ue)))    
 
    
           

