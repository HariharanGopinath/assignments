function [y]= rangekutta1euler(dt,y,t,f_x)                                                 
k = zeros(1,1);
A3=[0];
B3=[1];
C3=[0];
t=20;
        
for i=1:(length(t)-1)  
    
    k(1) = f_x(t(i));
          

    y(i+1) = y(i) + B3*k.'*dt; 
end