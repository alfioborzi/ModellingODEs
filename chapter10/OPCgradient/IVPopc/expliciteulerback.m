% Explicit Euler scheme backwards 
function [p] = expliciteulerback (N,h,yT,alpha,beta,y,u,B,a)
p(N+1)=-alpha*(y(N+1)-yT);

for n=(N+1):-1:2
t=(n-1)*h;
p(n-1)=p(n)+h*(dfdy(t,y(n),u(n),B,a)*p(n)-beta*(y(n)-yd(t)));
end
end
