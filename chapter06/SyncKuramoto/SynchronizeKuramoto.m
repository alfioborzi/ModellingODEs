
 % Numerical simulation for the Kuramoto model
%
% theta_i'=Omega_i + K/N sum_j=1^N sin(theta_j-theta_i)
%
clear all; close all;

ltype = {'b-','r--','m-.','k:'};             % for Plot
 
N=1000;   % number of oscillators
 
K=1.0;    % Coupling coefficient
 
h=0.1;
 
nSteps=500;
 
t=0:h:h*nSteps;
 
% initialize the theta_i with a uniform random distribution in 2pi
theta=zeros(N,nSteps);
theta(:,1)=2*pi*rand(N,1);

% initialize the phase speeds 
% omega=randn(N,1); K_c = 1.6
omega=(1-2*rand(N,1));  % K_c = 1.27

% Integration of the Kuramoto model with 4-th order RK method.
for j=1:nSteps
 
k1=kuramoto(theta(:,j),K,N,omega);
 
k2=kuramoto(theta(:,j)+0.5*h*k1,K,N,omega);
 
k3=kuramoto(theta(:,j)+0.5*h*k2,K,N,omega);           
 
k4=kuramoto(theta(:,j)+h*k3,K,N,omega);

% time step 
theta(:, j+1)=theta(:,j)+(h/6)*(k1+2*k2+2*k3+k4);
 

% plot the theta on the unit circle
x=cos(theta(:,j));
y=sin(theta(:,j));
 
s=linspace(0,2*pi,100);
cx=cos(s);
cy=sin(s);

plot(x,y,'o',cx,cy)
 
axis([-1 1 -1 1])
axis square
drawnow

end

% plots

% initial theta distribution 
figure(1)
% plot the theta on the unit circle
j = 1; 
x=cos(theta(:,j));
y=sin(theta(:,j));
 
s=linspace(0,2*pi,100);
cx=cos(s);
cy=sin(s);

plot(x,y,'o',cx,cy)
 
axis([-1 1 -1 1])
axis square

print('-depsc2', 'KuramotoInitTh01.eps','-b0'); 
print('-dpdf', 'KuramotoInitTh01.pdf','-b0');


% final theta distribution 
figure(2)
% plot the theta on the unit circle
j = nSteps; 
x=cos(theta(:,j));
y=sin(theta(:,j));
 
s=linspace(0,2*pi,100);
cx=cos(s);
cy=sin(s);

plot(x,y,'o',cx,cy)
 
axis([-1 1 -1 1])
axis square

print('-depsc2', 'KuramotoFinalTh01.eps','-b0'); 
print('-dpdf', 'KuramotoFinalTh01.pdf','-b0');



% order parameter 
p=zeros(size(t)); 

for j=1:nSteps+1
% order parameter 
realpart = 1/N*sum(arrayfun(@(y) cos(y),theta(:,j))); 
imagpart = 1/N*sum(arrayfun(@(y) sin(y),theta(:,j))); 
p(j) = sqrt((realpart)^2+(imagpart)^2);
end


figure(3)
plot(t,p,ltype{1},'Linewidth',1); 


axis([0 50 0 1])

print('-depsc2', 'orderPar01.eps','-b0'); 
print('-dpdf', 'orderPar01.pdf','-b0');
