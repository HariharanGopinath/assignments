clc;clear
close all

%% 2)a)
% Simulation parameters:
t_Final  = 2;
dt = 0.1;
dt_exact = 0.1;
dtRK2   = 0.1;
%From the question, assuming the lambda value
lambda=-2;
%initial value is one
x_0 = 1;


%itterations 
N=t_Final/dt;

%time period to calculate the solution
timeexact_sol=0:dt_exact:t_Final;
timeRK_sol=0:dt:t_Final;


% Exact solution
x_Exact = x_0*exp(lambda*timeexact_sol);
plot(timeexact_sol,x_Exact ,'marker','.','markersize',20)
hold on
legend('Exact solution')

% Explicit Euler Method:
% To find the euler true solution,we need to find the test function
syms x 
%the test function is 
x_dot=(lambda*x);

%converting that into a matlab function with x
f=matlabFunction(x_dot,'Vars',{x});

%butcher arrays for different integrators in the euler methods are. 
c= 0 ;
b= 1 ;
a= 0 ;        
%RK1
%epmty vectors for itterating
x_RK1 = zeros(N,1);
%the initial value is zero
x_RK1(1)=1;
for j = 2:N+1
    %evaluating x by euler method. 
    x_RK1(j) = x_RK1(j-1) + dt*f(x_RK1(j-1));
    
end
plot(timeRK_sol,x_RK1,'marker','.','markersize',10)
hold on;

c_2= [0 0.5] ;
b_2= [0 1] ;
a_2= [0,0;1/2 0] ;

%epmty vectors for itterating
x_RK2 = zeros(N,1);
%the initial value is zero
x_RK2(1)=1;

for j = 2:N+1
    %evaluating x by RK2 method.
    K1=f(x_RK2(j-1));
    K2=f(x_RK2(j-1)+(dt*a_2(2,1)*K1));
    x_RK2(j)=x_RK2(j-1)+(dt*b_2(1)*K1)+(dt*b_2(2)*K2);
end
plot(timeRK_sol,x_RK2,'marker','.','markersize',30)
hold on


%For the Rk4, the matrix values are
c_4= [0 0.5 0.5 1] ;
b_4= [1/6 1/3 1/3 1/6] ;
a_4= [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0 ] ; 

%epmty vectors for iterating
x_RK4 = zeros(N,1);
%the initial value is zero
x_RK4(1)=1;

for j = 2:N+1
    %evaluating x by RK4 method.
    K1=f(x_RK4(j-1));
    K2=f(x_RK4(j-1)+(dt*a_4(2,1)*K1));
    K3=f(x_RK4(j-1)+(dt*a_4(3,2)*K2));
    K4=f(x_RK4(j-1)+(dt*a_4(4,3)*K3));
    x_RK4(j)=x_RK4(j-1)+(dt*b_4(1)*K1)+(dt*b_4(2)*K2)+(dt*b_4(3)*K3)+(dt*b_4(4)*K4);
end
plot(timeRK_sol,x_RK4,'marker','.','markersize',40)
hold on
legend('Exact solution','RK1','RK2','RK4')
xlabel('Time')
ylabel('x')
title('Testing the Exact solution with Euler, RK2, RK4')
hold off
    
%% 2)b)
%Integration for various time steps:
N_table = floor(logspace(1,3,50));

for index = 1:length(N_table)
    
    N = N_table(index);
    dt(index) = t_Final / N_table(index);

    c_2= [0 0.5] ;
    b_2= [0 1] ;
    a_2= [0,0;1/2 0] ;
    %For the Rk4, the matrix values are
    c_4= [0 0.5 0.5 1] ;
    b_4= [1/6 1/3 1/3 1/6] ;
    a_4= [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0 ] ; 
    x_RK1 = zeros(N,1);
    x_RK2 = zeros(N,1);
    x_RK4 = zeros(N,1);
    %initial conditions
    x_RK1(1)=1;
    x_RK2(1)=1;
    x_RK4(1)=1;
 
    for j = 2:N+1
        
        %evaluating x by Rk1(euler method)
        x_RK1(j) = x_RK1(j-1) + dt(index)*f(x_RK1(j-1));
             
        %evaluating x by RK2 method.
        K1=f(x_RK2(j-1));
        K2=f(x_RK2(j-1)+(dt(index)*a_2(2,1)*K1));
        x_RK2(j)=x_RK2(j-1)+(b_2(1)*dt(index)*K1+b_2(2)*dt(index)*K2);
        
        %evaluating x by RK4 method
        
        K1=f(x_RK4(j-1));
        K2=f(x_RK4(j-1)+(dt(index)*a_4(2,1)*K1));
        K3=f(x_RK4(j-1)+(dt(index)*a_4(3,2)*K2));
        K4=f(x_RK4(j-1)+(dt(index)*a_4(4,3)*K3));
        x_RK4(j)=x_RK4(j-1)+dt(index)*((b_4(1)*K1)+(b_4(2)*K2)+(b_4(3)*K3)+(b_4(4)*K4));
        
        
    end
    
    % Compute the global error
    err(1,index) = norm(x_RK1(end)-x_Exact(end),inf);
    err(2,index) = norm(x_RK2(end)-x_Exact(end),inf);
    err(3,index) = norm(x_RK4(end)-x_Exact(end),inf);
end
figure(2); 
loglog(dt,err,'marker','.','markersize',15,'linestyle','none')
grid on
set(gca,'XDir','reverse')
legend('RK1', 'RK2','RK4')
xlabel('dt')
ylabel('Global error')
title('Global error vs dt')
    
%% c) For what value of l < 0 will the different schemes become unstable
lambda=(-2:-0.5:-40);
dt= 0.1;
tFinal=2;
N=t_Final/dt;
err_lambda=zeros(3,length(lambda));
unstable_RK1 = 0;
unstable_RK2 = 0;
unstable_RK4 = 0;
for index = 1:length(lambda)
    %test function symbolic
    xdot=lambda(index)*x;
    f=matlabFunction(xdot,'Vars',{x});
    %divide the total number of steps into equal siszed steps.
    x_RK1 = zeros(N,1);
    x_RK2 = zeros(N,1);
    x_RK4 = zeros(N,1);
    x_RK1(1)=1;
    x_RK2(1)=1;
    x_RK4(1)=1;
    for j = 2:N+1
        %evaluating x by Rk1(euler method)
        x_RK1(j) = x_RK1(j-1) + dt*f(x_RK1(j-1));
        
        %evaluating x by RK2 method.
        K1=f(x_RK2(j-1));
        K2=f(x_RK2(j-1)+(dt*a_2(2,1)*K1));
        x_RK2(j)=x_RK2(j-1)+(b_2(1)*dt*K1+b_2(2)*dt*K2);
        
        %evaluating x by RK4 method
        
        K1=f(x_RK4(j-1));
        K2=f(x_RK4(j-1)+(dt*a_4(2,1)*K1));
        K3=f(x_RK4(j-1)+(dt*a_4(3,2)*K2));
        K4=f(x_RK4(j-1)+(dt*a_4(4,3)*K3));
        x_RK4(j)=x_RK4(j-1)+(b_4(1)*dt*K1)+(b_4(2)*dt*K2)+(b_4(3)*dt*K3)+(b_4(4)*dt*K4);
        
        
    end
    
    % Compute the global error
    err_lambda(1,index) = norm(x_RK1(end)-x_Exact(end),inf);
    err_lambda(2,index) = norm(x_RK2(end)-x_Exact(end),inf);
    err_lambda(3,index) = norm(x_RK4(end)-x_Exact(end),inf);

    %checking error to determine the unstable condition.  
    %this formula can be found from the lecture notes for the stability
    %condition
    
   
    %stability check  for the RK1 
    if abs((lambda(index)*0.1)/1) <= 1
        unstable_RK1 = lambda(index);
    end
    
   %stability check  for the RK2  
   
    if  abs((lambda(index)*0.1) + ((lambda(index)*0.1)^2)/2) <= 1
        unstable_RK2 = lambda(index);
    end
    
    %stability check  for the RK4
    if abs(((lambda(index)*0.1) + ((lambda(index)*0.1)^2)/2 + (((lambda(index)*0.1)^3)/6) + (((lambda(index)*0.1)^4)/24))) <=1
        unstable_RK4 = lambda(index);
    end
end
disp(lambda)
figure(3)
plot(lambda,err_lambda)
legend('error  for RK1','error  for RK2','error for RK4')
xlabel('lambda')
ylabel('error')
title('error for different lambda')












    