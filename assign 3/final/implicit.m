clc; 
clear all; 
close all;
%% 4)a)
syms x dt real

%Defining the parameters and converting the symbol to a anonymous function

tfinal=1; 
lambda=-20; 
x_im = sym('x','real');
x_dot=lambda*x_im;
alp = 1;

x_dot_f=matlabFunction(x_dot,'Vars',{x});

tolerance=1e-5;
s=2;
n=1;


%Plothelp=linspace(0, tfinal, tfinal/dt)';

%Butcher
C = [(1/2)-(sqrt(3)/6);
    (1/2)+(sqrt(3)/6)] ;
A = [1/4 (1/4)-(sqrt(3)/6);
    (1/4)+(sqrt(3)/6) 1/4] ;
B = [1/2 1/2];  


K=sym('K',[s,n],'real');


%Here compute the necessary functions for the Newton method:

r= sym('r', [s,n],'real');
for i=1:s
    %evealuating the k 
    r(i,:)=x_dot_f(x_im + dt*(A(i,1)*K(1,:)+A(i,2)*K(2,:)))-K(i,:); 
end



dr=jacobian(r,K);
matlabFunction(r, dr, 'File','x_dot_r_dr', 'Vars',{K,x_im,dt})
clear dt r dr

dt=0.1;
N=tfinal/dt;
Plothelp=0: tfinal: N;


%Initial value as 1 and 1 in the first two column
xk = [1;1];
x_im = [xk;zeros(N-2,1)];

%Here we compute the Newton iteration:
for i=1:N-1
    %guess the k value 
    K = xk;
    
    while true
        [r,dr] = x_dot_r_dr(K,x_im(i),dt);
        norm(r)
        
        %taking the dr on the right hand side and it inverse of that
        dK = -dr\r;
        %update
        K = K + alp*dK;
        
        
        %break the loop
        if norm(r) < tolerance
            break
        end
      
            
    end
    K = reshape(K,[s,n]);
    x_im(i+1,:) = x_im(i,:) + dt*(B(1)*K(1)+B(2)*K(2));
end
Plothelp=linspace(0, tfinal, tfinal/dt)';


%% 4)b)
%Here we compute the explicit RK4 from the previous assignment:


%Defining the different symbols here:
syms x

%Defining the parameters and converting the symbol to a anonymous function
dt=0.1;
tfinal=1;
% lambda=-32; 
N=tfinal/dt;
x_dot=lambda*x;
fun=matlabFunction(x_dot,'Vars',{x});

%Compute the true solution:
x_0=1;
t_exact=(0:dt:tfinal)';
x_exact=x_0*exp(lambda*t_exact);
x_exact(N+1)=[];
%xEX(N+1)=[];


%Here compute the Euler:

c_4= [0 0.5 0.5 1] ;
b_4= [1/6 1/3 1/3 1/6] ;
a_4= [0 0 0 0; 
      1/2 0 0 0;
      0 1/2 0 0;
      0 0 1 0] ;
xRK4=ones(N,1);
%xRK4 = zeros(N,1);
%the initial value is zero
xRK4(1)=1;
for i=1:N-1
    k1=fun(xRK4(i));
    k2=fun(xRK4(i)+(dt*a_4(2,1))*k1);
    k3=fun(xRK4(i)+(dt*a_4(3,2))*k2);
    k4=fun(xRK4(i)+(dt*a_4(4,3))*k3);

    xRK4(i+1)=xRK4(i)+dt*((b_4(1)*k1)+(b_4(2)*k2)+(b_4(3)*k3)+(b_4(4)*k4));
    
end

figure(1) 
plot(Plothelp,x_im,'K',Plothelp,xRK4,'r',Plothelp,x_exact,'g');
legend('implicit','explicit','exact')
title('IRK4 VS ERK4 VS exact');
hold on;

% figure(2) 
% plot(Plothelp,xRK4,'r');
% title('ERK4');
% hold on
% 
% figure(3) 
% plot(Plothelp,x_exact,'r');
% title('exact');
% hold on


%errRK4=sqrt(mean((xRK4 - x_exact).^2));
%errx_RK4IM=sqrt(mean((x_RK4 - x_exact).^2));
%Diff=sqrt(mean((x_RK4 - xRK4).^2));%(Imp)

clear dt tfinal N

%% 4)c)
%symbolic variables
syms t x y
u=5;

%the nonlinear dynamics:
%xdot_Ydot=[y; u*(1-x^2)*y-x];

%Now creating a matlab function to convert the symbolic funtion in a matlab function
%xdot_Ydot_f = matlabFunction(xdot_Ydot,'Vars',{t,[x;y]}); 

%t=1;
% let assume the no of evaluation N 
%N=tf/dt;
%tfinal=25;
%dt=0.01;
%N=tfinal/dt;
%tRK_sol=zeros(N,1);
%tRK_sol(1)=1;
%tRK_sol=tRK_sol(2):dt:tfinal;


clc;
clear all;
close all;
%Defining the different symbols here:
syms dt t real

tolerance=1e-8;
s=2;
n=2;

%Defining the parameters and converting the symbol to a anonymous function

K = sym('K',[n,s],'real');

x = sym('x',[n,1]);    


xdot_Ydot=[x(2);
            5*(1-x(1)^2)*x(2) - x(1)];
        
xdot_Ydot_f = matlabFunction(xdot_Ydot,'Vars',{x}); 

%f1=y;
%f2=((u*((1-x_im.^2)*y))-x_im);
%y_f=matlabFunction(xdot_Ydot,'Vars',{x_im,y});

%x_f=matlabFunction(f2);
%Plothelp=linspace(0, tf ,N)';



% %Here we compute the Implicit RK4:
% ButchTIM=[1/2 - sqrt(3)/6, 1/4, 1/4-sqrt(3)/6;
%     1/2 + sqrt(3)/6, 1/4 + sqrt(3)/6, 1/4;
%     0, 1/2, 1/2];
% aRK4IM=ButchTIM(1:2,2:3);  
% bRK4IM=ButchTIM(3,2:3);
% cRK4IM=ButchTIM(1:2,1);

%y_van=zeros(N,1);
%y_van(1)=0

% Kim_van=sym('Kim_van',[n,s],'real');
% r1=sym('r1',[1,s*n],'real')';
% Kim_van=reshape(Kim_van,1,s*n);

%Butcher
C=[1/2-sqrt(3)/6;
   1/2+sqrt(3)/6];
B=[1/2;1/2];
A=[1/4, 1/4-sqrt(3)/6;
   1/4+sqrt(3)/6, 1/4];

% K1=sym('K',[1,s*n],'real');
%K1=sym('K',[1,s],'real');
%K2=sym('K',[1,s],'real');


%Here compute the necessary functions for the Newton method:

%r= sym('r', [1,s],'real');
%r= sym('r', [s*n,1],'real');

for i=1:s
    r(:,i)=xdot_Ydot_f(x + dt*(A(i,1)*K(:,1)+A(i,2)*K(:,2)))-K(:,i); 
end

%reshaping the values
r_change = reshape(r,[n*s,1]);
K_change = reshape(K,[n*s,1]);
% the jacobian can be found by
dr = jacobian(r_change,K_change);

func = matlabFunction(r_change, dr,'Vars',{t,K,x,dt});

clear dt t x r dr 
% for i=1:s
%     r(i,:)=xdot_Ydot_f(x_van + dt*(A(i,1)*K2(1,:)+A(i,2)*K2(2,:)), y_van+ dt*(A(i,1)*K2(1,:)+A(i,2)*K2(2,:)))-K2(i,:); 
% end


%for j = 1:s
%r1(:,j)=xdot_Ydot_f(x+dt*(aRK4IM(j,1)*Kim_van(1)+aRK4IM(j,2)*Kim_van(3)),y+dt*(aRK4IM(j,1)*Kim_van(2)+aRK4IM(j,2)*Kim_van(4)))-Kim_van(n*j-(n-1))-;
%end


%butcher array for the rk 4 from the last problem. 

tfinal=25;
dt=0.01;
N = tfinal/dt;

alpha = 1;
t=1;

c_4= [0 0.5 0.5 1] ;
b_4= [1/6 1/3 1/3 1/6] ;
a_4= [0 0 0 0; 
      1/2 0 0 0;
      0 1/2 0 0;
      0 0 1 0  ] ;
  
x_ini = [1;0];
x = [x_ini,zeros(n,N-1)];
K = [x_ini,x_ini];

for i = 1:N
    count = 0;
    alpha = 1;
    while true
        K=reshape(K,[n,s]);
        [r,dr] = func(t,K,x(:,i),dt);
      
        dk = -dr\r;

        K = reshape(K,[s*n,1]);

        K = K + alpha*dk;
        
        %To break out of the loop
        if norm(r) < tolerance
            break
        end
    end
    K = reshape(K,[n,s]);
    x(:,i+1) = x(:,i) + dt*(B(1)*K(:,1)+B(2)*K(:,2));

end


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
    
    K1 = xdot_Ydot_f([x_RK4(j-1);y_RK4(j-1)]);
    K2 = xdot_Ydot_f([x_RK4(j-1)+dt*0.5*K1(1);y_RK4(j-1)+dt*0.5*K1(2)]);
    K3 = xdot_Ydot_f([x_RK4(j-1)+dt*0.5*K2(1);y_RK4(j-1)+dt*0.5*K2(2)]);
    K4 = xdot_Ydot_f([x_RK4(j-1)+dt*1*K3(1);y_RK4(j-1)+dt*1*K3(2)]);
   %to find the XRK4 solution
    x_RK4(j) = x_RK4(j-1) + ( dt*(1/6) * K1(1) + dt*(1/3) * K2(1) + dt*(1/3)*K3(1) + dt*(1/6)*K4(1) );
    %to find the YRK4 solution
    y_RK4(j) = y_RK4(j-1) + ( dt*(1/6) * K1(2) + dt*(1/3) * K2(2) + dt*(1/3)*K3(2) + dt*(1/6)*K4(2) );
end

tRK_sol=linspace(0, tfinal, (tfinal/dt)+1)';
figure(2)
plot(tRK_sol,x(1,:))
hold on
plot(tRK_sol,x(2,:))
hold on

plot(tRK_sol,x_RK4)
hold on
plot(tRK_sol,y_RK4)
hold on
% figure(3)

legend('x1irk', 'x2irk','x_RK4','y_RK4')








