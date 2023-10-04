function t=genera_Poisson(lamb)
%
%  Simple way to generate Poisson density
%  pdf(t)=lamb*exp(-lamb*t)
%
%
unif = 0;
while (unif==0)
    unif = rand;
end
t = -log(unif)/lamb; 