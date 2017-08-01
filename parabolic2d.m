%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%             Parabolic 2D Problem                 %%%
%%%   u_t = u_xx + u_yy, (x,y) = (0,1)x(0,1), t > 0  %%% 
%%%    u(x,y,0) = sin(pi*x) sin(pi*y)                %%%
%%%     u = 0 on boundary                            %%%
%%% The exact soln: exp(-2*pi^2*t)sin(pi*x)sin(pi*y) %%%
%%%          Course: MATH F422, Dr. P. Dhanumjaya    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

a=0; b=1;  
c=0;   d=1;  
tfinal = 1.0;

n = 10;
   
   m = n;
   h = (b-a)/n;        
   dt=0.001;   
   x=a:h:b; 
   y=c:h:d;
   mu = dt/(2*h^2);

%-- Initial condition:

   for i=1:m+1,
      for j=1:m+1,
         u1(i,j) = sin(pi*x(i))*sin(pi*y(j)); 
      end
   end

%---------- loop for time t --------------------------------------

Tf = fix(tfinal/dt);

t = 0;

for k=1:Tf

for i = 1:m+1
    u2(i,1) = 0;
    u2(i,n+1) = 0;
end

A = sparse(m-1,m-1); b1=zeros(m-1,1);

for j = 2:n,                             % Look for fixed y(j) 

 for i=2:m,
      b1(i-1) = mu*(u1(i,j-1) -2*u1(i,j) + u1(i,j+1))+ ...
	 u1(i,j);
 end
      
 A(1,1) = (1+2*mu);
 A(1,2) = -mu;
 A(2,1) = -mu;
 for i = 2:m-2
     A(i,i+1) = -mu;
     A(i+1,i) = -mu;
     A(i,i) = (1+2*mu);
 end
 A(m-1,m-2) = -mu;
 A(m-1,m-1) = (1+2*mu);
 
     ut = A\b1;                          % Solve the diagonal matrix.
     for i=1:m-1,
	u2(i+1,j) = ut(i);
     end

 end                                    % Finish x-sweep.



%-------------- loop in y -direction --------------------------------

for i = 1:m+1
    u1(1,i) = 0;
    u1(m+1,i) = 0;
end

B = sparse(m-1,m-1); b2=zeros(m-1,1);

for i = 2:m,

   for j=2:n,
      b2(j-1) = mu*(u2(i-1,j) -2*u2(i,j) + u2(i+1,j))+ u2(i,j);
 end
      
 B(1,1) = (1+2*mu);
 B(1,2) = -mu;
 B(2,1) = -mu;
 
 for j = 2:n-2
     B(j,j+1) = -mu;
     B(j+1,j) = -mu;
     B(j,j) = (1+2*mu);
 end
 B(n-1,n-2) = -mu;
 B(n-1,n-1) = (1+2*mu);
 

     ut = B\b2;
     for j=1:n-1,
        u1(i,j+1) = ut(j);
     end

 end                             

 t = t + dt;
     
end       


  for i=1:m+1,
    for j=1:n+1,
       ue(i,j) = exp(-2*pi^2*t)*sin(pi*x(i))*sin(pi*y(j)); 
  end
  end
  
  error = max(max(abs(u1-ue)))    
 
  figure(1); 
  surf(x,y,u1);    
  title('Approximate solution')
  figure(2); 
  surf(x,y,ue)          
  title('Exact solution')
  