function [r,dr] = x_dot_r_dr(in1,x,dt)
%X_DOT_R_DR
%    [R,DR] = X_DOT_R_DR(IN1,X,DT)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    21-Oct-2020 10:59:42

K1 = in1(1,:);
K2 = in1(2,:);
r = [-K1-x.*2.0e1-dt.*(K1.*(1.0./4.0)-K2.*3.867513459481287e-2).*2.0e1;-K2-x.*2.0e1-dt.*(K1.*5.386751345948129e-1+K2.*(1.0./4.0)).*2.0e1];
if nargout > 1
    t2 = dt.*-5.0-1.0;
    dr = reshape([t2,dt.*(-1.077350269189626e1),dt.*7.735026918962573e-1,t2],[2,2]);
end