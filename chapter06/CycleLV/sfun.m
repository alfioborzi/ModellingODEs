function [sys] = sfun(t,x,mu)
% sfun
% van der Pool
%sys = [x(2), -x(1)+mu*x(2)*(1-x(1)^2)];

% Lotka - Volterra
 sys = [x(1)*(1-x(2)), mu*x(2)*(x(1)-1)];

% Limit cycle 
% sys = [x(2) + mu*(1-x(1)^2-x(2)^2)*x(1), -x(1) + mu*(1-x(1)^2-x(2)^2)*x(2)];

end