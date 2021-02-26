clc; clear all; close all;

%Defining the different symbols here:
syms x real

%Defining the parameters and converting the symbol to a anonymous function
dt=0.1;
tf=2; 
lamd=-2; 
N=tf/dt;
f1=lamd*x;
f=matlabFunction(f1,'Vars',{'x'});
toldo=1e-5;
Stages=2;
n=1;
Plothelp=linspace(0, tf, tf/dt)';

%Here we compute the Implicit RK4:
ButchTIM=[1/2 - sqrt(3)/6, 1/4, 1/4-sqrt(3)/6;1/2 + sqrt(3)/6, 1/4 + sqrt(3)/6, 1/4;0, 1/2, 1/2];
aRK4IM=ButchTIM(1:2,2:3);  
bRK4IM=ButchTIM(3,2:3);
cRK4IM=ButchTIM(1:2,1);
KIM_int=zeros(2,1)';
xRK4IM=ones(N,1);
KIM=sym('KIM',[Stages,n],'real')';

%Here compute the necessary functions for the Newton method:
for j = 1:Stages
    r1(j:j)=KIM(j)-f(x+dt*(aRK4IM(j,1)*KIM(1)+aRK4IM(j,2)*KIM(2)));
end
dr=jacobian(r1',KIM);
Fun_r1=matlabFunction(r1','Vars',{'KIM1','KIM2','x'});
Fun_r2=matlabFunction(dr\r1','Vars',{'KIM1','KIM2','x'});

%Here we compute the Newton iteration:
for i=1:N-1
    HELP_VAR=2;
    
    while HELP_VAR > 1
        KIM_int=KIM_int-Fun_r2(KIM_int(1),KIM_int(2),xRK4IM(i));
        
        if norm(Fun_r1(KIM_int(1),KIM_int(2),xRK4IM(i)))< toldo
            HELP_VAR=1;
        end
            
    end
    
    xRK4IM(i+1)=xRK4IM(i)+dt*((bRK4IM(1)*KIM_int(1))+(bRK4IM(2)*KIM_int(2)));
end

figure(1) 
plot(Plothelp,xRK4IM,'k');title('IRK4');
hold on


%Here we compute the explicit RK4 from the previous assignment:


%Defining the different symbols here:
syms x

%Defining the parameters and converting the symbol to a anonymous function
dt=0.1;
tf=2;
lamd=-2; 
N=tf/dt;
f=lamd*x;
ht=matlabFunction(f);

%Compute the true solution:
tEX=(0:dt:tf)';
xEX=exp(lamd*tEX);
xEX(N+1)=[];


%Here compute the Euler:
ButchT1=[0 0;0 1];
xE1=ones(N,1);
ButchT4 = [0 0 0 0 0; 0.5 0.5 0 0 0 ; 0.5 0 0.5 0 0; 1 0 0 1 0; 0 1/6 1/3 1/3 1/6];
aRK4 = ButchT4(1:4, 2:5);
bRK4 = ButchT4(5, 2:5);
cRK4 = ButchT4(1:4,1);
xRK4=ones(N,1);

for i=1:N-1
    Krk1=ht(xRK4(i));
    Krk2=ht(xRK4(i)+(dt*(sum(aRK4(2,:))*Krk1)));
    Krk3=ht(xRK4(i)+(dt*(sum(aRK4(3,:))*Krk2)));
    Krk4=ht(xRK4(i)+(dt*(sum(aRK4(4,:))*Krk3)));

    xRK4(i+1)=xRK4(i)+dt*((bRK4(1)*Krk1)+(bRK4(2)*Krk2)+(bRK4(3)*Krk3)+(bRK4(4)*Krk4));
    
end


%Plott them here: 

figure(2) 
plot(Plothelp,xRK4,'r');title('ERK4');

%Here we compute the true solution and the RMSE for both IRK4 and ERK4.
tEX=(0:dt:tf)';
xEX=exp(lamd*tEX);
xEX(N+1)=[];
ErrorRK4=sqrt(mean((xRK4 - xEX).^2));
ErrorxRK4IM=sqrt(mean((xRK4IM - xEX).^2));
Diff=sqrt(mean((xRK4IM - xRK4).^2));
%% Gauss-Legendre IRK (c)
clc;clear all; close all;
%Defining the different symbols here:
syms x y real

%Defining the parameters and converting the symbol to a anonymous function
u=5;
tf=25;
dt=0.01;
N=tf/dt;
f1=y;
f2=((u*((1-x.^2)*y))-x);
ht1=matlabFunction(f1,'Vars',{x,y});
ht2=matlabFunction(f2);
Plothelp=linspace(0, tf,N)';
toldo=1e-5;
Stages=2;
n=2;

%Here we compute the Implicit RK4:
ButchTIM=[1/2 - sqrt(3)/6, 1/4, 1/4-sqrt(3)/6;1/2 + sqrt(3)/6, 1/4 + sqrt(3)/6, 1/4;0, 1/2, 1/2];
aRK4IM=ButchTIM(1:2,2:3);  
bRK4IM=ButchTIM(3,2:3);
cRK4IM=ButchTIM(1:2,1);
KIM_int=zeros(n*Stages,n*Stages)';
xRK4IM=ones(N,1);
yRK4IM=ones(N,1);
KIM=sym('KIM',[n,Stages],'real');
r1=sym('r1',[1,Stages*n],'real')';
KIM=reshape(KIM,1,Stages*n);

%Here compute the necessary functions for the Newton method, and fill in r1:
for j = 1:Stages
    r1(n*j-(n-1):n*j-(n-1))=KIM(n*j-(n-1))-ht1(x+dt*(aRK4IM(j,1)*KIM(1)+aRK4IM(j,2)*KIM(3)),y+dt*(aRK4IM(j,1)*KIM(2)+aRK4IM(j,2)*KIM(4)));
end

for j = 1:Stages
     r1(n*j-(n-1)+1:n*j-(n-1)+1)=KIM(n*j-(n-1)+1)-ht2(x+dt*(aRK4IM(j,1)*KIM(1)+aRK4IM(j,2)*KIM(3)),y+dt*(aRK4IM(j,1)*KIM(2)+aRK4IM(j,2)*KIM(4)));
end
dr=jacobian(r1',KIM);
Fun_r1=matlabFunction(r1,'Vars',{'KIM1_1','KIM1_2','KIM2_1','KIM2_2','y','x'});
Fun_r2=matlabFunction(dr\r1,'Vars',{'KIM1_1','KIM1_2','KIM2_1','KIM2_2','y','x'});

%Here we compute the Newton iteration:
for i=1:N-1
    HELP_VAR=2;
    while HELP_VAR > 1
        KIM_int=KIM_int-Fun_r2(KIM_int(1),KIM_int(3),KIM_int(2),KIM_int(4),yRK4IM(i),xRK4IM(i));
        if norm(Fun_r1(KIM_int(1),KIM_int(3),KIM_int(2),KIM_int(4),yRK4IM(i),xRK4IM(i)))< toldo
            HELP_VAR=1;
        end     
    end
    
    yRK4IM(i+1)=yRK4IM(i)+dt*((bRK4IM(1)*KIM_int(2))+(bRK4IM(2)*KIM_int(4)));
    xRK4IM(i+1)=xRK4IM(i)+dt*((bRK4IM(1)*KIM_int(1))+(bRK4IM(2)*KIM_int(3)));
end
xRK4IM_y=[xRK4IM yRK4IM];

figure(1)
plot(Plothelp,xRK4IM_y);title('Van Der Pol for IRK4');
xlim([0 25])
ylim([-9 9])

%Here we compute the explicit RK4 from the previous assignment:
ButchT4 = [0 0 0 0 0; 0.5 0.5 0 0 0 ; 0.5 0 0.5 0 0; 1 0 0 1 0; 0 1/6 1/3 1/3 1/6];
aRK4 = ButchT4(1:4, 2:5);
bRK4 = ButchT4(5, 2:5);
cRK4 = ButchT4(1:4,1);
xRK4=ones(N,1);
y=ones(N,1);

  
  
for i=1:N-1
    Krk1=ht1(xRK4(i),y(i));
    Krk1dot=ht2(xRK4(i),y(i));
    
    Krk2=ht1(xRK4(i)+(dt*(sum(aRK4(2,:))*Krk1)),y(i)+(dt*(sum(aRK4(2,:))*Krk1dot)));
    Krk2dot=ht2(xRK4(i)+(dt*(sum(aRK4(2,:))*Krk1)),y(i)+(dt*(sum(aRK4(2,:))*Krk1dot)));
    
    
    Krk3=ht1(xRK4(i)+(dt*(sum(aRK4(3,:))*Krk2)),y(i)+(dt*(sum(aRK4(3,:))*Krk2dot)));
    Krk3dot=ht2(xRK4(i)+(dt*(sum(aRK4(3,:))*Krk2)),y(i)+(dt*(sum(aRK4(3,:))*Krk2dot)));
    
    Krk4=ht1(xRK4(i)+(dt*(sum(aRK4(4,:))*Krk3)),y(i)+(dt*(sum(aRK4(4,:))*Krk3dot)));
    Krk4dot=ht2(xRK4(i)+(dt*(sum(aRK4(4,:))*Krk3)),y(i)+(dt*(sum(aRK4(4,:))*Krk3dot)));

    xRK4(i+1)=xRK4(i)+dt*((bRK4(1)*Krk1)+(bRK4(2)*Krk2)+(bRK4(3)*Krk3)+(bRK4(4)*Krk4));
    y(i+1)=y(i)+dt*((bRK4(1)*Krk1dot)+(bRK4(2)*Krk2dot)+(bRK4(3)*Krk3dot)+(bRK4(4)*Krk4dot));
end
xRK4_y=[xRK4 y];

%Plott them here: 
figure(2)
plot(Plothelp,xRK4_y);title('Van Der Pol for ERK4');
xlim([0 25])
ylim([-9 9])

%Here we compute the RMSE between IRK4 and ERK4.
Diff=sqrt(mean((xRK4IM_y - xRK4_y).^2));