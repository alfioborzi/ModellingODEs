% opt_cont_elliptic solves a problem of the form
%   min J(x,u) := 1/2*||y-z||_L2^2 + nu/2*||u||_L2^2
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
function [u,y] = opc_string(z, u0, nu, f, tol,maxiter)
%
iter = 0;
u = u0;
n = length(u);
h = 1/(n-1);
x = linspace(0,1,n).';
f_h = f(x);
f_h = f_h(2:n-1);
z_h = z(x);
z_in = z_h(2:n-1);

% the -d^2 / d x^2 operator 
A_h = 1/h^2*(2*eye(n-2)-diag(ones(1,n-3),1)-diag(ones(1,n-3),-1));

% compute gradient in u0
gradJ_u = compgrad (u);
% evaluate functional J
[J_u,~] = compj(u);

fprintf('\n  k |  alpha | ||gradJ(u_k)|| | J\n');

% Method of steepest descent with Armijo rule
while sqrt(h)*norm(gradJ_u) > tol && iter < maxiter
  iter = iter + 1;
  alpha = linesearch(u,J_u, gradJ_u,h);
  u = u - alpha*gradJ_u;
  gradJ_u = compgrad (u);
  [J_u,x_n] = compj(u);
  fprintf('%3d | %.4f |  %6e |  %6e \n', iter, alpha, sqrt(h)*norm(gradJ_u),J_u);
end  

% Function to compute grad J(u)
% grad J(u) = nu*u - p
function [gradJ_u] = compgrad (u)
  u_in = u(2:n-1);
  y = A_h\(u_in+f_h);
  p = A_h\(-y+z_in);
  y = [0;y;0];
  p = [0;p;0];
  gradJ_u = nu*u - p;
end

% Function to compute J(u)
function [J_u,y] = compj(u)
  u_in = u(2:n-1);
  y = A_h\(u_in+f_h);
  y = [0;y;0];
  J_u = 0.5*h*norm(y-z_h)^2+nu/2*h*norm(u)^2;
end

% Linesearch with Armijo rule
% find alpha such that
% J(u-alpha*gradJ(u)) <= J(u)-alpha*d*||gradJ(u)||^2
function [alpha] = linesearch(u, J_u, gradJ_u,h)
  alpha = 1.0;
  d = 1.0e-2;
  norm2_grad = h*norm(gradJ_u);
  for k = 1:10
    J_new = compj(u-alpha*gradJ_u);   
    
    if J_new < ( J_u - alpha*d*norm2_grad ) 
        
      return;
    end
    alpha = alpha/2.;
  end
end

end