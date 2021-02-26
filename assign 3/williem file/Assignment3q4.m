%%a)----------------------------------------------------------------------
clear all 
close all
clc
%change this depending on wich function your going to run. 
%% switch for functions and conditions. 
subAssignment= 'c';
switch subAssignment
    case 'b'  
    %constants
    n=1;
    %Defining the test function:
    landa=-22;
    x  = sym('x','real');
    syms t u real

    f=landa*x;

    matlabFunction(f, 'file', 'f','vars',{x,u});
    clear u x t f
     
    case 'c'   
    %constants
    n = 2;
    %Defining VanDerPols:
    x  = sym('x',[n,1]);
    syms t u real

    f = [ x(2);
          u*(1-x(1)^2)*x(2) - x(1)];

    matlabFunction(f, 'file', 'f','vars',{x,u});
    clear x t f
end
%% Defining r(K,x,u): -----------------------------------------------------
%IRK method.
s=2;
% %symbolic vars. 
x = sym('x',[n,1]);
syms t dt u K real
K = sym('K',[n,s],'real');

%Butcher
C=[1/2-sqrt(3)/6;
   1/2+sqrt(3)/6];
B=[1/2;1/2];
A=[1/4, 1/4-sqrt(3)/6;
   1/4+sqrt(3)/6, 1/4];

%deploy K as a matrix. %now its hard coded do this in a loop. K(:,i) etc.
%for i in n*S. ADD u dependancy with c

%r=zeros(n,s);
for i=1:s
%Ks returns r=[n*s]
r(:,i)=f(x+dt*(A(i,1)*K(:,1)+A(i,2)*K(:,2)),u)-K(:,i); 
end

%K2
%r(:,2)=fv(x+dt*(A(2,1)*K(:,1)+A(2,2)*K(:,2)),u)-K(:,2);
%this reshaping causes errors all over the place.
%Reshape K and r to make the jacobian work
rr=reshape(r,[n*s,1]);
K_r=reshape(K,[n*s,1]);
dr = jacobian(rr,K_r);

%save function r
matlabFunction(rr,dr, 'file', 'rFileIRK4','vars',{t,K,x,dt,u});
clear xNext x t dt r dr
%% ------------------------------------------------------------------------

switch subAssignment
    case 'b'
    %order
    n=1;        
    %final time
    tFinal=1;
    %stepsize
    dt=0.01;
    %initial condition
    x0=[1];
    %guess for K?
    K=[x0,x0];
    %number of steps
    Nsteps=tFinal/dt;
    case 'c'
    %u
    u=5;
    %order
    n=2;
    %final time
    tFinal=25;
    %stepsize
    dt=0.1;
    %initial condition
    x0= [1;0];
    %guess for K?
    K=[x0,x0];  
    %number of steps
    Nsteps=tFinal/dt;
end
%assign t?
t=1;
%Create an empy vecotr for simulation
x = [x0,zeros(n,Nsteps-1)];
% Loop for the integrator with newton. 

for k = 1:Nsteps
    % Newton iteration
    iter = true;
    alpha = 0.5;
    niter = 0;
    while iter
         K=reshape(K,[n,s]); %reshape K back again, K has to go in r on this form. 
        [r,dr] = rFileIRK4(t,K,x(:,k),dt,u);
        
        dK=-dr\r; %delta K (newton step for deriving the sulotion for implicit equation in r with respect 2 K)
        K=reshape(K,[n*s,1]);
        K = K + alpha*dK;
        
        norm(r);
        if norm(r) < 1e-5
            iter = false;
        else
            niter = niter + 1;
        end
        
        if niter > 30
            disp('more then 30 loops something went wrong')
            iter=false;
        end
        
    end
    K=reshape(K,[n,s]);
    %for i=(1:s) %hardcoded for now, should be able to loop for n*s here
        x(:,k+1) = x(:,k)+dt*(B(1)*K(:,1)+B(2)*K(:,2));
end

%% ERK4 -------------------------------------------------------------------------------------
%renaming
stepsize=dt;
N=Nsteps;
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
    
    K1 = f([xRK4(j-1);yRK4(j-1)],u);
    K2 = f([xRK4(j-1)+stepsize*aRK4(2,1)*K1(1);yRK4(j-1)+stepsize*aRK4(2,1)*K1(2)],u);
    K3 = f([xRK4(j-1)+stepsize*aRK4(3,2)*K2(1);yRK4(j-1)+stepsize*aRK4(3,2)*K2(2)],u);
    K4 = f([xRK4(j-1)+stepsize*aRK4(4,3)*K3(1);yRK4(j-1)+stepsize*aRK4(4,3)*K3(2)],u);
    xRK4(j) = xRK4(j-1) + stepsize *( bRK4(1) * K1(1) + bRK4(2) * K2(1) + bRK4(3)*K3(1) + bRK4(4)*K4(1));
    yRK4(j) = yRK4(j-1) + stepsize *( bRK4(1) * K1(2) + bRK4(2) * K2(2) + bRK4(3)*K3(2) + bRK4(4)*K4(2));
end

%% -------------------------------------------------------------------------------------------------
%plotting for comparrison. 
plot(tPlot,x(1,:))
hold on
plot(tPlot,x(2,:))
hold on
plot(tPlot,xRK4)
hold on
plot(tPlot,yRK4)
hold off
legend('x1IRK4','x2IRK4','x1ERK4','x2ERK4')