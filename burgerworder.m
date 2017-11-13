%%non linear
%%u_t + u*u_x = u_xx
%%exact soln: ?


clear all;
close all;
clc
%exact solution for N=320
dt=0.0008;
tfin=1;
tini=0;
M=fix((tfin-tini)/dt);
tol=1e-8;
epsilon=3;
count=0;
N=320;
h=1/N;
F=zeros(N-1,1);
C=zeros(N-1,N-1);
U0=zeros(N-1,1);
U1=zeros(N-1,1);
F=zeros(N-1,1);



for j=1:N-1
    x(j)=h*j;
end
    for j=1:N-1
        U0(j)=sin(pi*x(j));
    end
    

mu=dt/(h^2);
lambda = dt/h ;
U=ones(N-1,1);
err=5;
t=0;

for k=1:M
      
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
        count = count +1;
    end
    U0=U1;
    t=t+dt;
    k;
    err = 5;
    %figure()
    %plot(x,U1,'o');
end
exact=zeros(N-1,1);
U1 = exact;

N=10;
h=1/N;
dt=0.0008;
tfin=1;
tini=0;
M=fix((tfin-tini)/dt);
tol=1e-8;
epsilon=3; %viscosity
k1=0;%no of iterations
t=0;
count=0;
for w=1:5 %space convergence
h=1/N;
F=zeros(N-1,1);
C=zeros(N-1,N-1);
U0=zeros(N-1,1);
U1=zeros(N-1,1);
F=zeros(N-1,1);
U10=zeros(N-1,1);
U20=zeros(N-1,1);
U30=zeros(N-1,1);
U40=zeros(N-1,1);
U50=zeros(N-1,1);


for j=1:N-1
    x(j)=h*j;
end
    for j=1:N-1
        U0(j)=sin(pi*x(j));
    end
    

mu=dt/(h^2);
lambda = dt/h ;
U=ones(N-1,1);
err=5;
t=0;

for k=1:M
      
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
        count = count +1;
    end
    U0=U1;
    t=t+dt;
    k;
    err = 5;
    %figure()
    %plot(x,U1,'o');
end


  if w ==1
      for k=1:N-1
          U10(k) = exact(k*32);
      end
      error(w) = max(abs(U10-U1))
  end
  if w == 2
     for k=1:N-1
     U20(k) = exact(16*k);
     end
     error(w) = max(abs(U20-U1))
  
  elseif w == 3
    for k=1:N-1
     U30(k) = exact(8*k);
     end
     error(w) = max(abs(U30-U1))
  
  elseif w == 4
     for k=1:N-1
     U40(k) = exact(4*k);
     end
     error(w) = max(abs(U40-U1))
  
  elseif w == 5
     for k=1:N-1
     U50(k) =exact(2*k);
     end
     error(w) = max(abs(U50-U1))
  end
  
N =2*N;
end
for p=2:4
  order(p) = log(error(p)/error(p+1))/log(2);
end
order




    

