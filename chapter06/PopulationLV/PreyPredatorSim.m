% Simulation of the modified Lotka-Volterra model
clear all; close all;

t = []; 		     % storage for time
x = []; 		     % storage for population size

t0 = 0; 

x0 = [0.2;0.4];     % Initial size of the population

a = 1.0; 
b = 0.5; 
d = 0.02; 

T = 100;

% Numerical solution
[t,x] = ode45(@(t,x) preypredatorODE(t,x,a,b,d) ,[t0 T],x0);


ltype = {'b-','r--','m-.'}; 


figure(1)
plot(t,x(:,1),ltype{1},t,x(:,2),ltype{2},'LineWidth',2);

xlabel('Time','Fontsize',14,...
        'Interpreter','LaTex');
ylabel('Population size',...
		'Fontsize',14,'Interpreter','LaTex');
legend({'$u$ - Prey','$v$ - Predator'},'Fontsize',14,...
        'Interpreter','LaTex');
    
print('-depsc2', 'modLV01.eps','-b0'); 
print('-dpdf', 'modLV01.pdf','-b0');

% Phase portrait
[X,Y] = meshgrid(0.1:0.1:1.0,0.1:0.1:1.0);
X = X';
Y = Y';
M = length(X(:,1));
N = length(X(1,:));
U = zeros(M,N);
V = zeros(M,N);
l = 0;
for m = 1:M
	for n = 1:N
		xdot = preypredatorODE(0,[X(m,n);Y(m,n)],a,b,d);
		l = 50*sqrt(xdot(1)^2+xdot(2)^2);
		U(m,n) = xdot(1)/l; % scaling
		V(m,n) = xdot(2)/l;	% scaling
	end
end

figure(2)
plot(x(:,1),x(:,2),ltype{3},'LineWidth',2);
hold on
quiver(X,Y,U,V,'AutoScale','off');
axis([0 0.8 0 0.6])
xlabel('Prey population size','Fontsize',14,...
        'Interpreter','LaTex');
ylabel('Predator population size','Fontsize',14,...
        'Interpreter','LaTex');
    
print('-depsc2', 'modLVcycle01.eps','-b0'); 
print('-dpdf', 'modLVcycle01.pdf','-b0');
