%One gradient step 
function [utemp,ytemp,Jtemp] = optimizeeuler (uold,dJ,y,y0,yT,alpha,beta,nu,N,h,B,a)

s=1;
% Parameter
sigma=10^-4;

% Evaluate the Functional in the old u and y
Jold=objective(y,uold,yT,alpha,beta,nu,N,h);

% compute a new value for u 
utemp=uold-s*dJ;
% solve the forward equation for this new u
ytemp=expliciteuler(h,N,y0,utemp,B,a);
% evaluate the functional in y and u
Jtemp=objective(ytemp,utemp,yT,alpha,beta,nu,N,h);

% line search
counter=0;
while Jtemp > Jold - s*sigma*sqrt(h)*norm(dJ)^2 && counter<=30
  s=s/2;
  % compute a new value for u   
  utemp=uold-s*dJ;  
  % solve the forward equation for this new u
  ytemp=expliciteuler(h,N,y0,utemp,B,a);
  % evaluate the functional in y and u
  Jtemp=objective(ytemp,utemp,yT,alpha,beta,nu,N,h);
  counter=counter+1;
end


end