%
% Code generate Poisson random variable
% H. Paul Keeler
%
function N=funPoissonRecursive(lambda)
T=0; %initialize sum of exponential variables as zero
n=-1; %initialize counting variable as negative one
%run (recursive) exponential function step function  

[~,N]=funStepExp(lambda,T,n);

function [T,N]=funStepExp(nu,S,m)
if(S<1)
%run if sum of exponential variables is not high enough

%generate exponential random variable
E=(-log(rand(1)))/nu;
S=S+E; %update sum of exponential variables
m=m+1; %update nunber of exponential variables

%recursively call function again
[T,N]=funStepExp(nu,S,m);
else
T=S;
N=m;
end
end
end
