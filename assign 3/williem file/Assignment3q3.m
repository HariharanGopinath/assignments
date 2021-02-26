
%% a)----------------------------------------------------------------------
clc;clear;
close all
%increasing the tollerance of ode45 to 1e-8. 
options = odeset('AbsTol',1e-8,'RelTol',1e-8);

%final time
tFinal=25;

%% Van der Polszz--------------------------- 
%symbolic variables
syms t x y
%constant
mu = 5;
vdp = [y; mu*(1-x^2)*y-x];
%create a function from the symbolic expression. 
vdpf = matlabFunction(vdp,'Vars',{t,[x;y]}); 
clear t x y
%% Get the stepsize from ode45---------------------------------------------
%run ode45 on the intervall t=[0,tFinal] for initial values
% x(0)=1, y(0)=0, 
%set up for ode45; ode45(function(f in LN),interval,initalconditions)
[ts,ys] = ode45(vdpf,[0 tFinal],[1;0]);
%stepsize in ode45
stepsize_ode45=zeros(length(ts),1);
    for i=1:length(ts)-1
    stepsize_ode45(i,1)=ts(i+1)-ts(i);
    end
%% b)------------------------------------------------------------------------

%RK4.
t=1; 
%stepsize=0.15;
N=1000;
%N=tFinal/stepsize;
stepsize=tFinal/N;

tPlot=zeros(N,1);
xRK4 = zeros(N,1);
yRK4 = zeros(N,1);
xRK4(1)=1;
yRK4(1)=0;
%butcher array
cRK4= [0 0.5 0.5 1] ;
bRK4= [1/6 1/3 1/3 1/6] ;
aRK4= [0 0 0 0; 
      1/2 0 0 0;
      0 1/2 0 0;
      0 0 1 0  ] ;
%itterations
for j = 2:N+1
    tPlot(j)=tPlot(j-1)+stepsize;
    
    K1 = vdpf(t,[xRK4(j-1);yRK4(j-1)]);
    K2 = vdpf(t,[xRK4(j-1)+stepsize*aRK4(2,1)*K1(1);yRK4(j-1)+stepsize*aRK4(2,1)*K1(2)]);
    K3 = vdpf(t,[xRK4(j-1)+stepsize*aRK4(3,2)*K2(1);yRK4(j-1)+stepsize*aRK4(3,2)*K2(2)]);
    K4 = vdpf(t,[xRK4(j-1)+stepsize*aRK4(4,3)*K3(1);yRK4(j-1)+stepsize*aRK4(4,3)*K3(2)]);
    xRK4(j) = xRK4(j-1) + stepsize *( bRK4(1) * K1(1) + bRK4(2) * K2(1) + bRK4(3)*K3(1) + bRK4(4)*K4(1));
    yRK4(j) = yRK4(j-1) + stepsize *( bRK4(1) * K1(2) + bRK4(2) * K2(2) + bRK4(3)*K3(2) + bRK4(4)*K4(2));
end


%Error in RK4 compared to ODE45 toll 1e-8
%what stepsize result in err <1e-8.
%cant get the error lower than this? error(0.01)=error(0.001) ??
err=zeros(1,2);
err(:,1) = norm(xRK4(end)-ys(end,1));
err(:,2)  = norm(yRK4(end)-ys(end,2));
% print the results in the command window:

fprintf('errx = %.6f \nerry = %.6f \n',err(:,1),err(:,2));

%RK4 compared to ODE45
figure(1)
plot(tPlot,xRK4)
hold on
plot(tPlot,yRK4)
hold on
plot(ts,ys(:,1))
hold on
plot(ts,ys(:,2))
hold off
legend('xRK4','yRK4','xODE45','yODE45')
title('VDP RK4 vs ODE45')
xlabel('x')
ylabel('time')

%stepsize RK4 VS ODE45
figure(2)
plot(ts,stepsize_ode45)
hold on
plot(ts,stepsize*ones(size(ts)))
hold off
legend('stepsize RK4','stepsize ODE45')
title('stepsize RK4 vs stepsize ODE45')
xlabel('time')
ylabel('stepsize')


