function [r_change,dr] = func(t,in2,in3,dt)
%FUNC
%    [R_CHANGE,DR] = FUNC(T,IN2,IN3,DT)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    22-Oct-2020 09:49:14

K1_1 = in2(1);
K1_2 = in2(3);
K2_1 = in2(2);
K2_2 = in2(4);
x1 = in3(1,:);
x2 = in3(2,:);
t2 = K2_1.*(1.0./4.0);
t3 = K2_2.*3.867513459481287e-2;
t4 = t2-t3;
t5 = dt.*t4;
t6 = K1_1.*(1.0./4.0);
t7 = K1_2.*3.867513459481287e-2;
t8 = t6-t7;
t19 = dt.*t8;
t9 = t19+x1;
t10 = K2_1.*5.386751345948129e-1;
t11 = K2_2.*(1.0./4.0);
t12 = t10+t11;
t13 = dt.*t12;
t14 = K1_1.*5.386751345948129e-1;
t15 = K1_2.*(1.0./4.0);
t16 = t14+t15;
t25 = dt.*t16;
t17 = t25+x1;
t18 = dt.*(1.0./4.0);
t20 = t5+x2;
t21 = t9.^2;
t22 = t21.*5.0;
t23 = t22-5.0;
t24 = dt.*5.386751345948129e-1;
t26 = t13+x2;
t27 = t17.^2;
t28 = t27.*5.0;
t29 = t28-5.0;
r_change = [-K1_1+t5+x2;-K2_1-x1-dt.*t8-t20.*t23;-K1_2+t13+x2;-K2_2-x1-dt.*t16-t26.*t29];
if nargout > 1
    dr = reshape([-1.0,-t18-dt.*t9.*t20.*(5.0./2.0),0.0,-t24-dt.*t17.*t26.*5.386751345948129,t18,dt.*t23.*(-1.0./4.0)-1.0,t24,dt.*t29.*(-5.386751345948129e-1),0.0,dt.*3.867513459481287e-2+dt.*t9.*t20.*3.867513459481287e-1,-1.0,-t18-dt.*t17.*t26.*(5.0./2.0),dt.*(-3.867513459481287e-2),dt.*t23.*3.867513459481287e-2,t18,dt.*t29.*(-1.0./4.0)-1.0],[4,4]);
end
