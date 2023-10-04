% Solve Cauchy problem
clear all;
close all;

ltype = {'b-*','m-.','k:','r--'};

% Solve ODE problem 
%
% y'   = e^y cos(x) 
% y(a) = ya
%

% odefun is f(x,y) in y'=f(x,y)
% odefunder is the derivative of odefun w.r.t. y
odefun = @(x,y) exp(y)*cos(x);
odefunder = @(x,y) exp(y)*cos(x);
ya = -1 ;

ee = exp(1);
sol = @(x) -log( ee - sin(x));

% interval [a,b]
a = 0; b = 2*pi;

% construct training data on a fine grid N > M 
Nx = 50; 
dx = (b-a)/Nx;

% x input values + bias - N points  
X = [a:dx:b]';
X = [X, ones(length(X),1)];

N = size(X,1);
 
% one input node + bias 
% one hidden layer with M nodes
% one output node 

M = 15;
alpha = 0.001;    % learning rate (eta)
maxepochs = 50000;   % number of epochs

    % activation function of hidden layer
    % l_sigm : logistic sigmoid
    % t_sigm : tanh 
    % r_sigm : ReLU
    
    act = @l_sigm;  

    % Work Units & Error
    WT = 0;
    ET = 0;
    WU = 0;
    
    disp('Learning iteration')
    
    % initializing weights from M grid 
    x = a + (b-a)*(0:(M))'/(M);
    dxM = (b-a)/M; 
    xi = [x, ones(length(x),1)];

    % All weights are learned
    % different initialization
    mcase = 1;  % fixed or random
    if mcase == 1 
    % with this setting one could train W2 only
    W1 = [2./(x(2:end)-x(1:(end-1))),(-1-2*x(1:(end-1))./(x(2:end)-x(1:(end-1))))];
    W2 = 0.0*[1:M]/M ;
    else    
    % random
    W1 = 2*rand(M, 2) -1;
    W2 = 2*rand(1, M) -1;
    end
    
    % learning
    
    % BP-SGD Algo
    for epoch = 1:maxepochs     
        
    % modified (simpler) BP
 %   [W1, W2, err] = BackpropODEmod(N, W1, W2, X, odefun, odefunder, 'alpha', alpha, 'a', a, 'ya', ya);

    % MLP-GCR Algo
    [W1, W2, err] = BackpropODEall(N, W1, W2, X, odefun, odefunder, 'alpha', alpha, 'a', a, 'ya', ya);

    % a WU=1 is a full BP step
    WU = WU + M/32 ;

    WT = [WT,WU];
    ET = [ET,err];
    end

    disp('Validation');
% exact solution on the grids N & M
        sN = sol(X(1:end,1));
        sM = sol(xi(1:end,1));
% NN solution on fine grid N
        zN = X(1:end,1); 
        for k = 1:N
          xr  = X(k,:)';
          v1 = W1*xr;
          y1 = act(0,v1);
          v  = W2*y1;
          y  = v; 
          zN(k) = TrialFunction(xr(1),y,'a',a,'ya',ya);
        end
% NN solution on M grid
        zM = xi(1:end,1); 
        for k = 1:M+1
          xM  = xi(k,:)';
          v1 = W1*xM;
          y1 = act(0,v1);
          v  = W2*y1;
          y  = v; 
          zM(k) = TrialFunction(xM(1),y,'a',a,'ya',ya);
        end        
        
% Num solution
        v1N = ode1(odefun,X(1:end,1),ya);
        v2N = ode2(odefun,X(1:end,1),ya);
        v1M = ode1(odefun,xi(1:end,1),ya);
        v2M = ode2(odefun,xi(1:end,1),ya);

           % Solutions on the M mesh 
            figure(2)
            plot(xi(1:end,1),zM,ltype{1},'Linewidth',2); hold on; 
            plot(xi(1:end,1),sM,ltype{2},'Linewidth',2); hold on; 
            plot(xi(1:end,1),v1M,ltype{3},'Linewidth',2);
            plot(xi(1:end,1),v2M,ltype{4},'Linewidth',2);
            hold off;
            legend('NN','Exact','1.order Num','2.order Num','Location','NorthEast');
             
 print('-depsc2', 'solutionsODEnn01.eps','-b0'); 
 print('-dpdf', 'solutionsODEnn01.pdf','-b0');
               
   
    errs = norm(zN-sN)*sqrt(dx);
    errn1 = norm(v1N-sN)*sqrt(dx);
    errn2 = norm(v2N-sN)*sqrt(dx);
    
    errsM = norm(zM(1:end,1)-sM(1:end,1))*sqrt(dxM);
    errn1M = norm(v1M(1:end,1)-sM(1:end,1))*sqrt(dxM);
    errn2M = norm(v2M(1:end,1)-sM(1:end,1))*sqrt(dxM);

    
    fprintf('N = %4i  errsN= %8.4f errn1N= %8.4f errn2N=%8.4f \n',N,errs,errn1,errn2)
    fprintf('M = %4i  errsM= %8.4f errn1M= %8.4f errn2M=%8.4f \n\n',M,errsM,errn1M,errn2M)

figure(3)
semilogy(WT,ET);
title('Loss function vs WU')

figure(4)
plot(W1(:,1),ltype{1},'Linewidth',2); hold on;
plot(W1(:,2),ltype{2},'Linewidth',2); hold on;
plot(W2(:),ltype{3},'Linewidth',2)
title('Weights')
legend('W1(:,1)','W1(:,2)','W2(:)','Location','NorthEast');

