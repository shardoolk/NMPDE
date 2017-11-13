%%non linear
%%u_t + u*u_x = u_xx
%%exact soln: ?


clear all;
close all;
clc
N=10;
h=1/N;
dt=0.008;
tfin=1;
tini=0;
M=fix((tfin-tini)/dt);
tol=1e-8;
epsilon=0; %viscosity
k1=0;%%no of iterations
t=0;
count=0;

A=zeros(N-1,N-1);
B=zeros(N-1,N-1);
U0=zeros(N-1,1);
U1=zeros(N-1,1);
F=zeros(N-1,1);


for j=1:N-1
    x(j)=h*j;
end
    for j=1:N-1
        U0(j)=sin(pi*x(j));
    end
    

mu=dt/(2*h^2);
lambda = dt/8*h ;
U=ones(N-1,1);
err=5;

for k=1:M
      
    while err>tol
         for j=1:N-1
            C(j,j)=1+2*mu*epsilon ;
        end
        for j=1:N-2
            C(j,j+1)=2*lambda*U(j+1) - epsilon*mu;
            C(j+1,j)=-2*lambda*U(j) - epsilon*mu;
        end
        F(1) = U(1)*(1 + 2*mu*epsilon) - U(2)*mu*epsilon - lambda*(U(1)^2) - (1-2*mu*epsilon)*U0(1) -  U0(2)*mu*epsilon - lambda*(U0(2)^2);
        F(N-1) = U(N-1)*(1 + 2*mu*epsilon) - U(N-2)*mu*epsilon - lambda*(U(N-1)^2) - (1-2*mu*epsilon)*U0(N-1) - U0(N-2)*mu*epsilon- lambda*(U0(N-2)^2);
        
        for j = 2:N-2
        F(j) = U(j)*(1 + 2*mu*epsilon) - U(j+1)*mu*epsilon -U(j-1)*mu*epsilon -lambda*(U(j+1)^2 - U(j-1)^2) - U0(j)*(1 - 2*mu*epsilon)- - U0(j+1)*mu*epsilon -U0(j-1)*mu*epsilon  -lambda*(U0(j+1)^2 - U0(j-1)^2);
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
    %figure()
    %plot(x,U1,'o');
end

plot(x,U1,'o');
            
