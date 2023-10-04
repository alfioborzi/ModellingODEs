% Derivative of the right hand side of the constraint ODE
function [dfdy] = dfdy (t,y,u,B,a)
dfdy=a*(1+t);
end
