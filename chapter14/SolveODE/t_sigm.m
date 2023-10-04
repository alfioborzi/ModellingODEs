function y = t_sigm(der,x)
% computes sigmoid function and its derivatives

switch(der)
    case 0
        y = tanh(x);
    case 1
        y = 1-t_sigm(0,x).*t_sigm(0,x);  
        y = diag(y);
    case 2
        y = -2.*t_sigm(0,x).*t_sigm(1,x);
         
end
end
