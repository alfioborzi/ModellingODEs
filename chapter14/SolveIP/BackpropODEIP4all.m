function [W1, W2] = BackpropODEIP4all(N, W1, W2, X, T)
  alpha = 0.1;
  
  e = zeros(4,1); 
 
  N1=N-1;
  
    x = X ;
    t = T ;
  
    v1 = W1*X(:);
    y1 = Sigmoid(v1);    
    v  = W2*y1;
    y  = v;       % these are Nj(W) , j=1,2,3,4
       
    % output of the model - derivatives
    for k = 1:N
    z(k,:)  = LotkaVolterraODE4all(x(k,:),y(1),y(2),y(3),y(4))  ;
    end
    
    % num derivative
    for k = 2:N1
    zp(k,:) = ((X(k+1,:)-X(k-1,:))/(t(k+1)-t(k-1))) ;
    end
    
    ea1 = 0.0;
    eb1 = 0.0;
    ec1 = 0.0;
    ed1 = 0.0;

    for k = 2:N1
    
    ea1 = ea1 + (  ( zp(k,1) - z(k,1) ) * x(k,1) ) ;
    eb1 = eb1 + ( - ( zp(k,1) - z(k,1) ) * (x(k,1)*x(k,2)) ) ;
    ec1 = ec1 + (  ( zp(k,2) - z(k,2) ) * (x(k,1)*x(k,2))  ) ;
    ed1 = ed1 + ( - ( zp(k,2) - z(k,2) ) * x(k,2) ) ;
    
    end
    
    e(1) = ea1 ;
    e(2) = eb1 ;
    e(3) = ec1 ;
    e(4) = ed1 ;
    
 % SGD & Backprop
 
    delta = e;

    e1     = W2'*delta;
    delta1 = y1.*(1-y1).*e1; 
    
    dW1 = alpha*delta1*X(:)';
    W1  = W1 + dW1;
    
    dW2 = alpha*delta*y1';    
    W2  = W2 + dW2;
  end
