function y = r_sigm(der,x)
% computes sigmoid function and its derivatives

switch(der)
    case 0
        y = max(0,x);
    case 1
        % an approx to ReLU is log(1 + e^x) 
        % its derivative is the logistic sigmoid
        y = l_sigm(0,x);
        y = diag(y);  
    case 2
        % with the approx above 
  %       y = l_sigm(0,x);
        y = l_sigm(0,x).*(1-l_sigm(0,x));  
        y = diag(y);        
end
end
