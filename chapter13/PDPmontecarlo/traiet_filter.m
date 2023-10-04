function [x,t]=traiet_filter(x0,t0,stato,dt,param)
%
%  [x,t]=traiet_filter(x0,t0,stato,dt,gamma,W)
%  Given initial config. (x0,t0), the dyn state 'stato' 
%  and dt, returns the position x at t
%  gamma and W parameters of the process
%
g = param.gam(stato);
w = param.W(stato);

dtn=linspace(0,dt,30);

decay = exp(-g*dtn);
x = w/g.*(1.-decay)+x0.*decay;
t=t0+dtn;

