
% opt_cont_elliptic solves a problem of the form
%   min J(y,u) := 1/2*||y - z||_L2^2 + nu/2*||u||_L2^2
%      s.t. - y'' = f + u  in (0,1)
%                    y = 0      on boundary
% by using the method of steepest descent with Armijo rule.
%
% Input:  z = target function
%        u0 = first guess of control u (discrete)
%         f = function of elliptic equation
%        nu = scalar nu
%
% Output: u = optimal control 
%         y = y(u)
%
% Given data:
% 

clear all; close all; 

ltype = {'b-','r--','m-.','k:'};             % for Plot

% control weight

nu = 1.0e-5;

% rhs and target
 f = @(x) 0.*x;
 z = @(x) sign(min(0,sin(3*pi*x)));
   
 
 % discretization
n = 200;                     % number of grid points
xvector = linspace(0,1,n).';

% max SD iterations and tol for the norm of red. gradient
maxiter = 10000;
tol = 10^-5;

%To test with an exact solution:
%  f   = @(x) pi^2*sin(pi*x)-1/nu*sin(2*pi*x);
%  z   = @(x) 4*pi^2*sin(2*pi*x)+sin(pi*x);
%  u_s = @(x) (1/nu)*sin(2*pi*x);
%  y_s = @(x) sin(pi*x);


u0 = zeros(n,1);            % starting iter with control = 0

% Optimization
[u, y] = opc_string(z, u0, nu, f, tol, maxiter);

% plot

figure(1)
z_h = z(xvector);
plot(xvector,z_h,ltype{1},'Linewidth',1);         % plot the target z

print('-depsc2', 'targetEll.eps','-b0'); 
print('-dpdf', 'targetEll.pdf','-b0');


figure(2)
plot(xvector,y,ltype{1},'Linewidth',1);         % plot the solution y

print('-depsc2', 'stateEll.eps','-b0'); 
print('-dpdf', 'stateEll.pdf','-b0');


figure(3)
plot(xvector,u,ltype{1},'Linewidth',1);         % plot the solution u

print('-depsc2', 'controlEll.eps','-b0'); 
print('-dpdf', 'controlEll.pdf','-b0');
