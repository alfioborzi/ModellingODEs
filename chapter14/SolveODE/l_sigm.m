function y = l_sigm(der,x)
% computes sigmoid function and its derivatives

switch(der)
    case 0
        y = 1./(1+exp(-x));
    case 1
        y = l_sigm(0,x).*(1-l_sigm(0,x));  
        y = diag(y);
    case 2
        y = l_sigm(1,x).*(1-2*l_sigm(0,x));
end
end
