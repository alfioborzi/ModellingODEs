% Reduced gradient of the optimal control problem
function [dJ] = grad (u,p,y,nu,B)
dJ=nu*u-B*p;
end