% The functional which is the objective function of the optimal Control
% Problem
%
% Input: Parameter alpha, beta, yT, nu 
%        Arguments y, u
%
function [J] = objective (y,u,yT,alpha,beta,nu,N,h)
sum1=0;
sum2=0;

for i=2:N
  t=(i-1)*h;
  sum1=sum1+(y(i)-yd(t))^2;
  sum2=sum2+u(i)^2;
end
sum2=sum2+u(N+1)^2+u(1)^2;
J=0.5*beta*(h*sum1)+0.5*alpha*(y(N+1)-yT)^2+0.5*nu*(h*sum2);
end