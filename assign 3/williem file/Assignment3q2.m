
%% a)-----------------------------------------------------------------------
%explicit euler and explicit RK2
clc;clear;
close all;clear all

%test function exact value, and plot in given interval.  
xzero=1;
stepsize= 0.1;
stepsizeExact=0.1;
landa=-2;
tFinal= 2;
tRK=0:stepsize:tFinal;
tExact = 0:stepsizeExact:tFinal;
xExact = xzero*exp(landa*tExact);
figure(1)
plot(tExact,xExact)
hold on
%symbolic variables
syms x xdot real
%test function symbolic
xdot=landa*x;

f=matlabFunction(xdot,'Vars',{x});

%itterations 
N=tFinal/stepsize;

%butcher arrays for different integrators. 
cRK1= [0] ;
bRK1= [1] ;
aRK1= [0] ;        
%RK1
%epmty vectors for itterating
xRK1 = zeros(N,1);
xRK1(1)=1;
for j = 2:N+1
    %evaluating x by euler method. 
    xRK1(j) = xRK1(j-1) + stepsize*f(xRK1(j-1));
end
%plotting the integrator
plot(tRK,xRK1)
legend('Rk1')
hold on
cRK2= [0 0.5] ;
bRK2= [0 1] ;
aRK2= [0,0;1/2 0] ;    
%RK2
%epmty vectors for itterating
xRK2 = zeros(N,1);
xRK2(1)=1;

for j = 2:N+1
    K1 = f(xRK2(j-1));
    K2 = f(xRK2(j-1)+aRK2(2,1)*stepsize*K1);
    xRK2(j) = xRK2(j-1) + stepsize * bRK2(1) * K1 + stepsize * bRK2(2) * K2 ;
                 
end
%plotting the integrator
plot(tRK,xRK2)
hold on
cRK4= [0 0.5 0.5 1] ;
bRK4= [1/6 1/3 1/3 1/6] ;
aRK4= [0 0 0 0; 
      1/2 0 0 0;
      0 1/2 0 0;
      0 0 1 0  ] ;   
%RK4
%epmty vectors for itterating
xRK4 = zeros(N,1);
xRK4(1)=1;

for j = 2:N+1
    K1 = f(xRK4(j-1));
    K2 = f(xRK4(j-1)+stepsize*aRK4(2,1)*K1);
    K3 = f(xRK4(j-1)+stepsize*aRK4(3,2)*K2);
    K4 = f(xRK4(j-1)+stepsize*aRK4(4,3)*K3);
    xRK4(j) = xRK4(j-1) + stepsize *( bRK4(1) * K1 + bRK4(2) * K2 + bRK4(3)*K3 + bRK4(4)*K4);
end
%plotting the integrator
plot(tRK,xRK4)
legend('Analytical sulotion','RK1','RK2','RK4')
xlabel('time')
ylabel('x')
title('test equation')
hold off

%% b)-------------------------------------------------------------------------------
%Error vs timestep. 

%Integration for various time steps:
%vector with different number of timesteps. 50 N from [10-1000]
N_table = floor(logspace(1,3,50));

%make a simulation for each N (number of timesteps)
for index = 1:length(N_table)
    
    N = N_table(index);
    
    %divide the total number of steps into equal siszed steps.
    stepsize(index) = tFinal / N_table(index);
    t = zeros(N,1);
    xRK1 = zeros(N,1);
    xRK2 = zeros(N,1);
    xRK4 = zeros(N,1);
    xRK1(1)=1;
    xRK2(1)=1;
    xRK4(1)=1;
    for j = 2:N+1
       
        % RK1
        xRK1(j) = xRK1(j-1) + stepsize(index)*f(xRK1(j-1));
        
        % RK2 
        
        K1 = f(xRK2(j-1));
        K2 = f(xRK2(j-1)+aRK2(2,1)*stepsize(index)*K1);
        xRK2(j) = xRK2(j-1) + stepsize(index) * bRK2(1) * K1 + stepsize(index) * bRK2(2) * K2 ;
         
        %RK4   
        K1 = f(xRK4(j-1));
        K2 = f(xRK4(j-1)+stepsize(index)*aRK4(2,1)*K1);
        K3 = f(xRK4(j-1)+stepsize(index)*aRK4(3,2)*K2);
        K4 = f(xRK4(j-1)+stepsize(index)*aRK4(4,3)*K3);
        xRK4(j) = xRK4(j-1) + stepsize(index) *( bRK4(1) * K1 + bRK4(2) * K2 + bRK4(3)*K3 + bRK4(4)*K4);

    end
    
    % Compute the global error for each stepsize(index) (global error is the local error at t=T=stepsize*N)
    err(1,index) = norm(xRK1(end)-xExact(end),inf);
    err(2,index) = norm(xRK2(end)-xExact(end),inf);
    err(3,index) = norm(xRK4(end)-xExact(end),inf);
end
%plot the log log of the different errors. 
figure(2); 
loglog(stepsize,err,'marker','.','markersize',15,'linestyle','none')
grid on
set(gca,'XDir','reverse')
legend('RK1', 'RK2','RK4')
xlabel('stepsize')
ylabel('error')
title('error vs stepsize logscale')

%% c)----------------------------------------------------------------------------------------
%stability check. 
landa=(-2:-0.5:-30);
stepsize= 0.1;
tFinal= 2;
N=tFinal/stepsize;
err=zeros(3,length(landa));
RK1_unstable = 0;
RK2_unstable = 0;
RK4_unstable = 0;
for index = 1:length(landa)
    %test function symbolic
    xdot=landa(index)*x;
    f=matlabFunction(xdot,'Vars',{x});
    %divide the total number of steps into equal siszed steps.
    xRK1 = zeros(N,1);
    xRK2 = zeros(N,1);
    xRK4 = zeros(N,1);
    xRK1(1)=1;
    xRK2(1)=1;
    xRK4(1)=1;
    for j = 2:N+1
       
        % RK1
        xRK1(j) = xRK1(j-1) + stepsize*f(xRK1(j-1));
        
        % RK2 
        
        K1 = f(xRK2(j-1));
        K2 = f(xRK2(j-1)+aRK2(2,1)*stepsize*K1);
        xRK2(j) = xRK2(j-1) + stepsize * bRK2(1) * K1 + stepsize * bRK2(2) * K2 ;
         
        %RK4   
        K1 = f(xRK4(j-1));
        K2 = f(xRK4(j-1)+stepsize*aRK4(2,1)*K1);
        K3 = f(xRK4(j-1)+stepsize*aRK4(3,2)*K2);
        K4 = f(xRK4(j-1)+stepsize*aRK4(4,3)*K3);
        xRK4(j) = xRK4(j-1) + stepsize *( bRK4(1) * K1 + bRK4(2) * K2 + bRK4(3)*K3 + bRK4(4)*K4);

    end
    
    % Compute the global error for each stepsize(index) (global error is the local error at t=T=stepsize*N)
    err(1,index) = norm(xRK1(end)-xExact(end),inf);
    err(2,index) = norm(xRK2(end)-xExact(end),inf);
    err(3,index) = norm(xRK4(end)-xExact(end),inf);
    %check the error to determine when the integrators become unstable. 

    if err(1,index) > 1
        if RK1_unstable == 0
        % RK 1 unstabel for 
        RK1_unstable=landa(index)
        end
    end
    if err(2,index) > 1
        if RK2_unstable == 0
        % RK 2 unstabel for 
        RK2_unstable=landa(index)
        end
    end
    if err(3,index) > 1
        if RK4_unstable == 0
        % RK 4 unstabel for 
        RK4_unstable=landa(index)
        end
    end
end

disp('running like a charm')
figure(3)
plot(landa,err)
legend('error RK1','error RK2','error RK4')
xlabel('landa')
ylabel('error')
title('error for different landa')