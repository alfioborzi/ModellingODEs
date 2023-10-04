function [W1, W2, E] = BackpropODEall(N, W1, W2, X, odefun, odefunder, varargin)
% learning SGD wise from data x such with loss function
% E = 0.5*(y_bar' - odefun(y_bar))^2
% with y_bar = ya + (x-a)*N(x,w);
% odefunder is the derivative of odefun w.r.t. y
    a = 0;
    ya = 0;
    alpha = 0.0001;
    
    act = @l_sigm;	% activation function of hidden layer, 1st argument is number of derivative
    
    for i = 1:2:length(varargin)
        switch varargin{i} 
            case 'a'
                a = varargin{i+1};
            case 'ya'
                ya = varargin{i+1};
            case 'alpha'
                alpha = varargin{i+1};
            case 'act'
                act = varargin{i+1};
        end
    end
    
    error = zeros(1,N);
    for k = 1:N
        x = X(k,:)';
        
        v1 = W1*x;
        y1 = act(0,v1);
        v2 = W2*y1;
        y2 = v2;                    % = N(x,w)
        
        y_bar = TrialFunction(x(1),y2,'a',a,'ya',ya);
        
        y2der = W2*act(1,v1)*W1(:,1);
        y_barder = DerivativeTrialFunc(x(1),y2,y2der,'a',a);
        
        e = -(y_barder - odefun(x(1),y_bar));
       
        dW2 = e*(1-odefunder(x(1),y_bar)*(x(1)-a))*y1' + e*(x(1)-a)*W1(:,1)'*act(1,v1);
        dW1 = e*((1-odefunder(x(1),y_bar)*(x(1)-a))*act(1,v1)*W2' + (x(1)-a)*W2'.*act(2,v1)*W1(:,1))*x' + e*(x(1)-a)*act(1,v1)*W2'*[1,0];
        
        W2 = W2 + alpha*dW2;
        W1 = W1 + alpha*dW1;
        
        error(k) = 0.5*e^2;
    end
    E = mean(error);
end
