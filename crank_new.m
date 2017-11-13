clear all
close all
clc
N=10;
h=1/N;
dt=0.001;
tfin= 0.2;
tini=0;
M=fix((tfin-tini)/dt);
tol=1e-8;
epsilon=0.1; %viscosity
k1=0;%%no of iterations
t=0;
count=0;

A=zeros(N-1,N-1);
C=zeros(N-1,N-1);
U0=zeros(N-1,1);
U1=zeros(N-1,1);
F=zeros(N-1,1);


for j=1:N-1
    x(j)=h*j;
        U0(j)=sin(pi*x(j));
    end
    

r=(epsilon*dt)/(2*h^2);
p = dt/4*h ;
U=ones(N-1,1);
err=5;
for k=1:M
      
    while err>tol
        F(1) = (-r + p*U0(1))*U(2) + (2*r+1 +p*U0(2))*U(1) - r*U0(2) - (1-2*r)*U0(1);
        for j = 2:N-2
            F(j) = (-r + p*U0(j))*U(j+1) + (2*r+1+ p*U0(j+1) - p*U0(j-1))*U(j) - (r + p*U0(j))*U(j-1) - r*U0(j+1) -(1-2*r)*U0(j) -r*U0(j-1);
        end
        F(N-1) = (2*r + 1 - p*U0(N-2))*U(N-1) -1*(r + p*U0(N-1))*U(N-2) - (1-2*r)*U0(N-1) -r*U0(N-2);
            C(1,1) = 1 + 2*r + p*(U0(2));
            for i = 2:N-2
                C(i,i) = 1 + 2*r + p*(U0(i+1) - U0(i-1));
            end 
            C(N-1,N-1) = 1 + 2*r + p*(U0(N-2));
            
            for j = 1:N-2
                C(j,j+1) = -r  + p*U0(j);
                C(j+1,j) = -r - p*U0(j+1);
            end 
            
        delta=C\(-F);
        U1=delta+U;
        U=U1;
        err=max(abs(delta));
        count = count +1
    end
    U0=U1;
    t=t+dt;
    k
    err = 5;
  % figure()
   %plot(x,U1,'o');
end
plot(x,U1,'o');
U1
    