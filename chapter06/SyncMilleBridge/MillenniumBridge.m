% Sync of Millennium Bridge

close all; clear all;

ltype = {'b-','r--','m-.','k:'};             % for Plot

% 
% y'' + b y' + omega0^2 y = g sum_j sin (theta_j)
%
% theta' = w_j + C A sin ( Psi - theta_j + alpha)
%

omega0 = 6.0;      % resonance frequency omega_0 = sqrt(k/m)
a = omega0^2;    
b = 0.1;              % drag b = 2 beta
g = 1.0*10^-3;        % coupling coefficient 

% Time setting
t0 = 0;          % time begin
t1 = 3600*2.0;        % time end  30 min. => t1 = 1800
% Time grid from t0 to T
T      = t1-t0;
Nsteps = 100000;
h      = T/Nsteps;

iTime = Nsteps+1;
Time  = t0:h:T;

% Initial conditions
x1 = 0.0;       % initial position of the bridge
v1 = 0.0;       % initial velocity of the bridge

% Pedestrian, number of oscillators
N=300;   
% initialize the phase theta_i with a uniform random distribution in 2pi
theta1=2*pi*rand(N,1);

% initialize the phase speeds of the oscillators
% omega=(1-2*randn(N,1)) + omega0;  
sigma = 1.0;
omega = sigma*randn(N,1) + omega0;

% phase lag
alpha = pi/2 ;
C = 20.;

% Initial conditions for the bridge and the pedestrian
wtheta1 = [x1; v1; theta1];

% Bridge and pedestrian dynamics  
[wtheta]=ode5(@(t,wtheta)KuramotoMFbridge(wtheta,omega,omega0,a,b,g,C,alpha),Time,wtheta1);



% figure(1)
% for j=1:30:Nsteps
% % plot the theta on the unit circle
% xth=cos(wtheta(j,3:end));
% yth=sin(wtheta(j,3:end));
%  
% s=linspace(0,2*pi,100);
% cx=cos(s);
% cy=sin(s);
% 
% plot(xth,yth,'o',cx,cy)
%  
% axis([-1 1 -1 1])
% axis square
% drawnow
% 
% end

% order parameter 
p=zeros(size(Time)); 
Nlast=Nsteps/10;     % 3 min

for j=1:Nsteps
% order parameter 
realpart = 1/N*sum(arrayfun(@(y) cos(y),wtheta(j,3:end))); 
imagpart = 1/N*sum(arrayfun(@(y) sin(y),wtheta(j,3:end))); 
p(j) = sqrt((realpart)^2+(imagpart)^2);
end

% plot order parameter of the pedestrian - Nlast steps
figure(1)
plot(Time(end-Nlast:end),p(end-Nlast:end),ltype{1},'Linewidth',1); 

print('-depsc2', 'orderPed01.eps','-b0'); 
print('-dpdf', 'orderPed01.pdf','-b0');


% plot lateral displacement of the bridge
figure(2)
plot(Time(end-Nlast:end),wtheta(end-Nlast:end,1),ltype{1},'Linewidth',0.1);
% hold on; 
% plot(Time,wtheta(:,2)); hold on; 

print('-depsc2', 'displaceB01.eps','-b0'); 
print('-dpdf', 'displaceB01.pdf','-b0');

