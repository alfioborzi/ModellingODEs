function f=kuramoto(x,K,N,omega)
 
f=omega+(K/N)*sum(sin(x*ones(1,N)-(ones(N,1)*x')))';
 
end