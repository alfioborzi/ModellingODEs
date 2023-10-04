% Stability test doe Euler-Maruyama Method
%
% 	dX = lambda*X dt + mu*X dW,   X(0) = Xzero,
%   
% with lambda = -3, mu = 2, and Xzero = 1.

rng(100,'v5normal')
T = 5; M = 10000; Xzero = 1;          
ltype = {'b-','r--','m-.'};             % for Plot

figure(1)
%%%%%%%%%%%% square-mean %%%%%%%%%%%%%
lambda = -3; mu = 2;              % Problem's parameter
for k = 1:3
    Dt = 2^(-k-1) ;                     
    N = T/Dt;
    Xms = zeros(1,N); Xtemp = Xzero*ones(M,1);
    for j = 1:N
           Winc = sqrt(Dt)*randn(M,1);  
           Xtemp = Xtemp + Dt*lambda*Xtemp + mu*Xtemp.*Winc;
           Xms(j) = mean(abs(Xtemp).^2);     
    end
 %   semilogy([0:Dt:T],[Xzero,Xms],ltype{k},'Linewidth',2), 
       plot([0:Dt:T],[Xzero,Xms],ltype{k},'Linewidth',2),
    hold on
end
legend({'$h = 0.25$','$h = 0.125$',...
	'$h = 0.0625$'},'Interpreter','Latex','FontSize',12)
%title('Stabilit{\"a}t im quadratischen Mittel','FontSize',16,...
%	'Interpreter','Latex')
%ylabel('$E\left(\left|X\right|^2\right)$','FontSize',14,...
%	'Interpreter','Latex')
axis([0,T,-0.25,1.75]), 
hold off

print('-depsc2', 'msStabEMSDE01.eps','-b0'); 
print('-dpdf', 'msStabEMSDE01.pdf','-b0');
