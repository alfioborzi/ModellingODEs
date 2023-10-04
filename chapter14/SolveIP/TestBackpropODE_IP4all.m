clear all;
close all;

ltype = {'k-','b-*','m-.','ro'};

% ODE problem - Lotka Volterra
%
% x' (1) = x(1) * ( a - b * x(2)) 
% x' (2) = x(2) * ( c * x(1) - d )

% Construct input data 

t = []; 		     % storage for time
x = []; 		     % storage for population size


t0 = 0; 
x0 = [0.2;0.3];     % Initial size of the population

% these parameters should be recovered from NN comp
a = 1.0;
b = 2.0; 
c = 1.0;
d = 0.3; 


% Time horizon, approx. one cycle 
N=100;
T0 = 13;
h  = (T0-t0)/N; 
t = linspace(t0,T0,N);

% Numerical solution

 x = ode5(@(t,x) LotkaVolterraODE4(x,a,b,c,d),t,x0);

% x input values - measuraments of population size with time
% take every second result of the num simulation

X = x(1:2:N,:);
T = t(1:2:N);

% number of training inputs 
Nt = size(X,1);
     
% two time Nt input nodes
% one hidden layer with K nodes
% four output nodes 

K  = 20;
Ni = 2*Nt;
W1 = eye(K, Ni) ;
W2 = eye(4, K) ;

for epoch = 1:50000          % train
  [W1 W2] = BackpropODEIP4all(Nt, W1, W2, X, T);
end

     
    v1 = W1*X(:);
    y1 = Sigmoid(v1);    
    v  = W2*y1;
    y  = v;     % these are Nj(W)    
    
    a1 = y(1);
    b1 = y(2);
    c1 = y(3);
    d1 = y(4);

disp('estimated parameter values ')

fprintf('a = %8.4f  b= %8.4f c= %8.4f  d=%8.4f \n',a1,b1,c1,d1)

% Figures 

% Numerical solution with estimated parameters 

 z = ode5(@(t,z) LotkaVolterraODE4(z,a1,b1,c1,d1),t,x0);
 
% Phase portrait
[Z,Y] = meshgrid(0.1:0.1:1.0,0.1:0.1:1.0);
Z = Z';
Y = Y';
M = length(Z(:,1));
N = length(Z(1,:));

% with a and b 
U = zeros(M,N);
V = zeros(M,N);
l = 0;
for m = 1:M
	for n = 1:N
		xdot = LotkaVolterraODE4([Z(m,n);Y(m,n)],a,b,c,d);
		l = 50*sqrt(xdot(1)^2+xdot(2)^2);
		U(m,n) = xdot(1)/l; % scaling
		V(m,n) = xdot(2)/l;	% scaling
	end
end



figure(1)
% compare x and z 
plot(x(:,1),x(:,2),ltype{1},'LineWidth',2); hold on; 
plot(z(:,1),z(:,2),ltype{4},'LineWidth',2); 
hold on
quiver(Z,Y,U,V,'AutoScale','off');
axis([0 0.8 0 1.0])
xlabel('Prey population size','Fontsize',14,...
        'Interpreter','LaTex');
ylabel('Predator population size','Fontsize',14,...
        'Interpreter','LaTex');
    
print('-depsc2', 'solLVcycleAll01.eps','-b0'); 
print('-dpdf', 'solLVcycleAll01.pdf','-b0');




