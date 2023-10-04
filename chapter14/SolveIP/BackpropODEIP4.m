function [W1, W2] = BackpropODEIP4(N, W1, W2, X, T)
  alpha = 0.1;
  
 e = zeros(4,1); 
  
  for k = 2:N-1
    x = X(k,:)';
    t = T(k);
    
    v1 = W1*x;
    y1 = Sigmoid(v1);    
    v  = W2*y1;
    y  = v;       % these are Nj(W) , j=1,2,3,4
    
    % output of the model - derivatives
    z  = LotkaVolterraODE4(x,y(1),y(2),y(3),y(4))  ;
    % num derivative
    zp = ((X(k+1,:)-X(k-1,:))/(T(k+1)-T(k-1)))' ;
    
    e(1) =   ( zp(1) - z(1) ) * x(1)  ;
    e(2) = - ( zp(1) - z(1) ) * (x(1)*x(2)) ;
    e(3) =   ( zp(2) - z(2) ) * (x(1)*x(2)) ;
    e(4) = - ( zp(2) - z(2) ) * x(2) ;
    
    delta = e;

    e1     = W2'*delta;
    delta1 = y1.*(1-y1).*e1; 
    
    dW1 = alpha*delta1*x';
    W1  = W1 + dW1;
    
    dW2 = alpha*delta*y1';    
    W2  = W2 + dW2;
  end
end