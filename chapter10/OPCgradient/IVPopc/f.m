% Righthand side of the ODE constraint
function [f] = f (y,t,u,B,a)
f=a*(1+t)*y+B*u;
end
