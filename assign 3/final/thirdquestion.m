clc;clear;
close all

% From, question, final time 
tf=25;

%to set the default options for ode45 
options = odeset('AbsTol',1e-8,'RelTol',1e-8);



%symbolic variables
syms t x y
u=5;

%the nonlinear dynamics:
xdot_Ydot=[y; (u*(1-x^2)*y)-x];

%Now creating a matlab function to convert the symbolic funtion in a matlab function
xdot_Ydot_f = matlabFunction(xdot_Ydot,'Vars',{t,[x;y]}); 

% now, the x and y observes the value from function
clear  x y
%% 3)a run ode45 on the above function
%defining the initial conditions
x_0y_0=[1;0];
%defining the time interval
tspan=[0 25];
%defining the ode function
[tf,yf]=ode45(xdot_Ydot_f,tspan,x_0y_0,options);

dt_ode45=zeros(length(tf),1);
    for i=1:length(tf)-1
    dt_ode45(i,1)=tf(i+1)-tf(i);
    end



%% 3)b Deploy your own RK4 scheme on the same dynamics, and compare the solutions.
t=1;
%let assume the N value 
N=1000;
tfinal=25;
dt=tfinal/N;
tRK_sol=zeros(N,1);
%tRK_sol=0:dt:tfinal;
%butcher array for the rk 4 from the last problem. 

c_4= [0 0.5 0.5 1] ;
b_4= [1/6 1/3 1/3 1/6] ;
a_4= [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];
  
x_RK4 = zeros(N,1);
y_RK4 = zeros(N,1);
%assigning the initial conditions
x_RK4(1)=1;
y_RK4(1)=0;

for j = 2:N+1
    
    

    %K1 = xdot_Ydot_f(t,[x_RK4(j-1);y_RK4(j-1)]);
    %K2 = xdot_Ydot_f(t,[x_RK4(j-1)+dt*a_4(2,1)*K1(1);y_RK4(j-1)+dt*a_4(2,1)*K1(2)]);
    %K3 = xdot_Ydot_f(t,[x_RK4(j-1)+dt*a_4(3,2)*K2(1);y_RK4(j-1)+dt*a_4(3,2)*K2(2)]);
    %K4 = xdot_Ydot_f(t,[x_RK4(j-1)+dt*a_4(4,3)*K3(1);y_RK4(j-1)+dt*a_4(4,3)*K3(2)]);
    %x_RK4(j) = x_RK4(j-1) + dt *( b_4(1) * K1(1) + b_4(2) * K2(1) + b_4(3)*K3(1) + b_4(4)*K4(1));
    %y_RK4(j) = y_RK4(j-1) + dt *( b_4(1) * K1(2) + b_4(2) * K2(2) + b_4(3)*K3(2) + b_4(4)*K4(2));
    %evaluating x by RK4 method.
    K1 = xdot_Ydot_f(t,[x_RK4(j-1);y_RK4(j-1)]);
    K2 = xdot_Ydot_f(t,[x_RK4(j-1)+dt*0.5*K1(1);y_RK4(j-1)+dt*0.5*K1(2)]);
    K3 = xdot_Ydot_f(t,[x_RK4(j-1)+dt*0.5*K2(1);y_RK4(j-1)+dt*0.5*K2(2)]);
    K4 = xdot_Ydot_f(t,[x_RK4(j-1)+dt*1*K3(1);y_RK4(j-1)+dt*1*K3(2)]);
   %to find the XRK4 solution
    x_RK4(j) = x_RK4(j-1) + ( dt*(1/6) * K1(1) + dt*(1/3) * K2(1) + dt*(1/3)*K3(1) + dt*(1/6)*K4(1) );
    %to find the YRK4 solution
    y_RK4(j) = y_RK4(j-1) + ( dt*(1/6) * K1(2) + dt*(1/3) * K2(2) + dt*(1/3)*K3(2) + dt*(1/6)*K4(2) );
end


%finiding the error solution comparing with the ode and rk solution
err_xdot_Ydot=zeros(1,2);
err_xdot_Ydot(:,1) = norm(x_RK4(end)-yf(end,1));

err_xdot_Ydot(:,2)  = norm(y_RK4(end)-yf(end,2));

disp(err_xdot_Ydot(:,1));
disp(err_xdot_Ydot(:,2))

tRK_sol=linspace(0, tfinal, (tfinal/dt)+1)';
% plot([1:length(tf)],tf)
%plotting the graphs 
%RK4 compared to ODE45
figure(1)
plot(tRK_sol,x_RK4)
hold on
plot(tRK_sol,y_RK4)
hold on
plot(tf,yf(:,1))
hold on
plot(tf,yf(:,2))
hold off
legend('xRK4','yRK4','xODE45','yODE45')
title('Van der pol RK4 solution vs ODE45 solution')
xlabel('x')
ylabel('time')

%Comparing the discrete time grid to the one of ode45.
figure(2)
plot(tf,dt_ode45)
hold on
plot(tf,dt*ones(size(tf)))
hold off
legend('dt ODE45', 'dt RK4')
title('Ploting the step size of RK4 solution vs step size of ODE45 solution')



