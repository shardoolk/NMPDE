%non linear
%u_t + u*u_x = u_xx
%exact soln: ?


clear all;
close all;
clc
N=80;
h=1/N;
dt=0.001;
tfin= 0.8;
tini=0;
M=fix((tfin-tini)/dt);
tol=1e-8;
epsilon=0.1 ;%viscosity
t=0;
count=0;

F=zeros(N-1,1);
C=zeros(N-1,N-1);
U0=zeros(N-1,1);
U1=zeros(N-1,1);
U2=zeros(N-1,1);
F=zeros(N-1,1);

for j=1:N-1
    x(j)=h*j;
    U0(j)=sin(pi*x(j));
end
    

mu=dt/(h^2);
lambda = dt/h ;
U=ones(N-1,1);
err=5;
%following is for getting value at t=dt 
while err>tol
         for j=1:N-1
            C(j,j)=1+2*mu*epsilon ;
        end
        for j=1:N-2
            C(j,j+1)=lambda*U(j+1)/2 - epsilon*mu;
            C(j+1,j)=-lambda*U(j)/2 - epsilon*mu;
        end
        F(1)=(1+2*epsilon*mu)*U(1) + (U(2)*lambda/4 -epsilon*mu)*U(2) - U0(1) ;
        F(N-1)=(1+2*epsilon*mu)*U(N-1) + (-U(N-2)*lambda/4 -epsilon*mu)*U(N-2) - U0(N-1) ;
        for j=2:N-2
            F(j)=(1+2*epsilon*mu)*U(j) + (U(j+1)*lambda/4 -epsilon*mu)*U(j+1) + (-U(j-1)*lambda/4 -epsilon*mu)*U(j-1) - U0(j);
        end
        
        
        delta=C\(-F);
        U1=delta+U;
        U=U1;
        err=max(abs(delta));
       
    end
t=dt;
mu=2*mu;
lambda = dt/(2*h) ;
U1 = zeros(N-1,1);
for k=1:M-1
    U1(1) = ( U(2)*(-U(2)*lambda + mu*epsilon) + U0(1)*(1 - mu*epsilon) )/(1 + mu*epsilon) ;  
    U1(N-1) = ( U(N-2)*(U(N-2)*lambda + mu*epsilon) + U0(N-1)*(1 - mu*epsilon) )/(1 + mu*epsilon) ;
    for p=2:N-2
        U1(p) =( U(p+1)*(-U(p+1)*lambda + mu*epsilon) +U(p-1)*(U(p-1)*lambda + mu*epsilon)+ U0(p)*(1 - mu*epsilon) )/(1 + mu*epsilon) ;
    end
    U0=U;
    U=U1;
    t=t+dt;
%   figure()
 %  plot(x,U1,'o');
end
  
    


for i=1:N-1
%   exact(i)=exp(-1)*sin(pi*x(i));
end

%dev=max(abs(exact'-U1));
%plot(x,U1,'o');
max(U1)


