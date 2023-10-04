% Solve the optimal control problem
function [y,p,u,normdJ,Jplot,counter] = optimalexpliciteuler (a,B,t0,T,N,alpha,beta,nu,y0,yT,epsilon,maxiter)

h=T/N;
u=zeros(1,N+1);


dJ=100000*ones(1,N+1);
counter=0;

% Forward
y=expliciteuler(h,N,y0,u,B,a);

Jplot(counter+1)=objective(y,u,yT,alpha,beta,nu,N,h);
while sqrt(h)*norm(dJ)>epsilon && counter<=maxiter

% Backward
p=expliciteulerback (N,h,yT,alpha,beta,y,u,B,a);
% Compute gradient
dJ=grad(u,p,y,nu,B);
% Walk one gradient step
[u,y,J]=optimizeeuler(u,dJ,y,y0,yT,alpha,beta,nu,N,h,B,a);


counter=counter+1;
normdJ(counter)=sqrt(h)*norm(dJ);
Jplot(counter+1)=J;

end
end
