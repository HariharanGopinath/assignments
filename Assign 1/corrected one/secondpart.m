clear all;
clc;

syms p1 p2 x1 y1 z1 x2 y2 z2 dx1 dy1 dz1 dx2 dy2 dz2 m1 m2 g l U ux uy uz Z;


% the generalized co ordinates,as per the question, will be
p1=[x1; y1; z1];
p2=[x2; y2; z2];
q=[p1;p2];
%the input force in this is same as last problem
U=[ux; uy; uz; 0 ;0; 0];


e=p1-p2;

C=0.5*(transpose(e)*e-(l^2));

p1_dot=[dx1; dy1; dz1];
p2_dot=[dx2; dy2; dz2];
dq=[p1_dot;p2_dot];
jacob_p1=jacobian(p1,q);
jacob_p2=jacobian(p2,q);
dp1=jacob_p1*dq;
dp2=jacob_p2*dq;

t1=0.5*m1*transpose(dp1)*dp1;
t2=0.5*m2*transpose(dp2)*dp2;
T=t1+t2;

% potential energies of both the helicopter and hanging masses are
v1=-m1*g*z1;
v2=-m2*g*z2;
V=v1+v2;

W1=m1*transpose(jacob_p1)*jacob_p1;
W2=m2*transpose(jacob_p2)*jacob_p2;
W=W1+W2;

Wdq=W*dq;
Jacobian_Wdq=jacobian(Wdq,q);
b2=U-(Jacobian_Wdq*dq)-jacobian(T,q).'-jacobian(V,q).'-Z*jacobian(C,q).';

%Solving 2a problem with same values from previous problem 
%here we need to define and specify the a and c 


a=jacobian(C,q);

%Now, we need find the function of C(dq,q,u)
jacobC_q=jacobian(C,q);
jacboC_dq=jacobian((jacobC_q*dq),q);

C_final=[(b2+Z*jacobian(C,q).'); -(jacboC_dq*dq)];

%now we need to solve the M defined in the question

M=[W transpose(jacobC_q); jacobC_q 0];

inverse_M=inv(M);

%Finally we need to find the answers with respect to dq_dot and Z

dq_dot=inverse_M*C_final;
A=latex(b2)
