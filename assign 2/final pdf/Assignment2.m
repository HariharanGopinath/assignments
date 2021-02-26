clc
clear all
load input
load output

% splitting the values
N = length(y);  % number of data
uestimation = u(1:N/2);
yestimation = y(1:N/2);
uvalidation = u(N/2+1:end);
yvalidation = y(N/2+1:end);

PHI = zeros(3,N/2);           
PHI(:,1) = [ 0 ; 0 ; uestimation(1)];
PHI(:,2) = [ yestimation(1) ; 0 ; uestimation(2) ];
for t=3:N/2 %considering half a value
    PHI(:,t) = [yestimation(t-1) ; yestimation(t-2) ; uestimation(t)] ;  
end


% first equation y(t)=-a1y(t-1)-a2y(t-2)+b0u(t)+e(t)

th = (PHI*PHI')\PHI*yestimation;

a11 = th(1);
a12 = th(2);
b10 = th(3);


%%


PHI = zeros(4,N/2);           
PHI(:,1) = [ 0 ; 0 ; uestimation(1) ; 0];
PHI(:,2) = [ yestimation(1) ; 0 ; uestimation(2) ; uestimation(1) ];
for t=3:N/2 %considering half a value
    PHI(:,t) = [yestimation(t-1) ; yestimation(t-2) ; uestimation(t) ; uestimation(t-1)] ;  
end



% second equation y(t)=-a1y(t-1)-a2y(t-2)+b0u(t)+b1u(t-1)+e(t)
th = (PHI*PHI')\PHI*yestimation;

a21 = th(1);
a22 = th(2);
b20 = th(3);
b21 = th(4);




%%



PHI = zeros(4,N/2);           
PHI(:,1) = [ 0 ; 0 ; 0 ; 0];
PHI(:,2) = [ yestimation(1) ; 0 ; 0 ; uestimation(1) ];
PHI(:,3) = [ yestimation(2) ; yestimation(1) ; 0 ; uestimation(2) ];

for t=4:N/2 %considering half a value
    PHI(:,t) = [yestimation(t-1) ; yestimation(t-2) ; yestimation(t-3) ; uestimation(t-1)] ;  
end



% Third equation y(t)=-a1y(t-1)-a2y(t-2)-a3y(t-3)+b1u(t-1)+e(t)
th = (PHI*PHI')\PHI*yestimation;

a31 = th(1)
a32 = th(2)
a33 = th(3)
b31 = th(4)


%%
ypredication1 = zeros(N/2,1);
ysimulation1 = zeros(N/2,1);

for t=3:N/2 %considering half a value
    ypredication1(t) = a11*yvalidation(t-1)+a12*yvalidation(t-2)+b10*uvalidation(t) ;
    ysimulation1(t) = a11*ysimulation1(t-1)+a12*ysimulation1(t-2)+b10*uvalidation(t) ;

end

predERROR1 = yvalidation-ypredication1;
predRMSE1  = rms(predERROR1)
simERROR1 = yvalidation-ysimulation1;
simRMSE1  = rms(simERROR1)


%%

ypred2 = zeros(N/2,1);
ysim2 = zeros(N/2,1);

for t = 4:N/2 %considering half a value
    ypred2(t)= a21*yvalidation(t-1)+a22*yvalidation(t-2)+b20*u(t)+b21*uvalidation(t-1) ;
    ysim2(t)= a21*ysim2(t-1)+a22*ysim2(t-2)+b20*u(t)+b21*uvalidation(t-1) ;
end

predERROR2 = yvalidation-ypred2;
predRMSE2  = rms(predERROR2)
simERROR2 = yvalidation-ysim2;
simRMSE2  = rms(simERROR2)

%%


ypred3 = zeros(N/2,1);
ysim3 = zeros(N/2,1);

for t = 4:N/2 %considering half a value
    ypred3(t)= a31*yvalidation(t-1)+a32*yvalidation(t-2)+a33*yvalidation(t-3)+b31*uvalidation(t-1);
    ysim3(t)= a31*ysim3(t-1)+a32*ysim3(t-2)+a33*ysim3(t-3)+b31*uvalidation(t-1);
end


predERROR3 = yvalidation-ypred3;
predRMSE3  = rms(predERROR3)
simERROR3 = yvalidation-ysim3;
simRMSE3  = rms(simERROR3)


%% Plotting the graph

figure (1)
subplot(3,1,1) %splitting the graph
plot(yvalidation)
hold on
plot(ypredication1)
legend('data points','prediction of model 1')
title('Output')
xlabel('Samples')
ylabel('output')
grid on

subplot(3,1,2)
plot(yvalidation)
hold on
plot(ypred2)
legend('data points','prediction of model 2')
title('Output')
xlabel('Samples')
ylabel('output')
grid on

subplot(3,1,3)
plot(yvalidation)
hold on
plot(ypred3)
legend('data points','prediction of model 3')
title('Output')
xlabel('Samples')
ylabel('output')
grid on 

figure (2)
subplot(3,1,1)
plot(yvalidation)
hold on
plot(ysimulation1)
legend('data points','simulation of modal 1')
title('Output')
xlabel('Samples')
ylabel('output')
grid on

subplot(3,1,2)
plot(yvalidation)
hold on
plot(ysim2)
legend('data points','simulation of modal 2')
title('Output')
xlabel('Samples')
ylabel('output')
grid on

subplot(3,1,3)
plot(yvalidation)
hold on
plot(ysim3)
legend('data points','simulation of modal 3')
title('Output')
xlabel('Samples')
ylabel('output')
grid on
