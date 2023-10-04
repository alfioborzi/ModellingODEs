% Explicit Euler Scheme 
function [y] = expliciteuler (h,N,y0,u,B,a)
y(1)=y0;

for n=1:N
t=(n-1)*h;
y(n+1)=y(n)+h*f(y(n),t,u(n),B,a);
end

end
