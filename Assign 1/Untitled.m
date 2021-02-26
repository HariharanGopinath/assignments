clear all;
clc;
syms x y z l q p1 p2 g t1 t2 v1 v2 m1 m2 dx dy dz dq p1dot p2dot teeta phi teeta_dot phi_dot ux uy uz;
% the generalized coordinates are
q=[p1; teeta; phi];
%expanding the p1(position of helicopter) with x y z coordinates
p1= [x; y; z];
% the newly formed generalized coordinates will be 
q=[x; y; z; teeta; phi];
% now the position of the hanging mass with respect to x y z will be
p2=[x+ l*cos(teeta)*sin(phi); y+ l*sin(phi)*sin(teeta); z-l*cos(phi)];

%the differentiation of p1 and p2
dp1=[dx; dy; dz];
dq=[dp1; teeta_dot; phi_dot];

dp2=jacobian(p2,q)*dq;

% kinetic energies of both the helicopter and hanging masses are
t1=0.5*m1*transpose(dp1)*dp1;
t2=0.5*m2*transpose(dp2)*dp2;
T=t1+t2;

% potential energies of both the helicopter and hanging masses are
v1=m1*g*z;
v2=m2*g*(z-l*cos(phi));
V=v1+v1;
 
%lagrange equation
L=T-V;

%Euler lagrange equations
% grad_Ldq=jacobian(L,dq).';
% diff_grad_Ldq=jacobian(grad_Ldq,q);
% grad_L=jacobian(L,q).';
% Q=diff_grad_Ldq-grad_L

%Now we need to calculate the W(q) for both the masses

jacob_p1=jacobian(p1,q);
jacob_p2=jacobian(p2,q);

W1=m1*transpose(jacob_p1)*jacob_p1;
W2=m2*transpose(jacob_p2)*jacob_p2;
W=W1+W2;

%From the question(M(q)?v = b(q;dq;u)), now we need to define the function
%of b with respect to q, dq and input parameter u.

%the input parameter U can be defined as
U=[ux; uy; uz; 0 ;0 ];

%the extended Euler-lagrange equation
W_dq= W*dq;
JacobianW_dq=jacobian(W_dq,q);
b=U-(JacobianW_dq*dq)-jacobian(T,q).'-jacobian(V,q).'







