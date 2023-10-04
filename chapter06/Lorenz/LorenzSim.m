%Simulation of the Lorenz model 
close all; clear all;       

% Typical Parameters
sigma = 10;
b = 8/3;
r = 28;

% The Loren system
f=@(t,x) lorenz(x,sigma,r,b);

t0 = 0; 		%initial time
T = 50; 		%final time
%Startwerte
x0 = [1;1;1];
y0 = [0.9999998;1.0000002;1.0000002];
 
t1 = []; 	    % storage 
t2 = [];
x1 = [];	    				 
x2 = [];


%Chaos dynamics : high-order accurate integration

[t1,x1]=  ode113(@(t,x) f(t,x),[t0 T],x0);
[t2,x2]=  ode113(@(t,x) f(t,x),[t0 T],y0);

% For axes
xMax = max(max(x1(:,1)),max(x2(:,1)));
yMax = max(max(x1(:,2)),max(x2(:,2)));
zMax = max(max(x1(:,3)),max(x2(:,3)));
xMin = min(min(x1(:,1)),min(x2(:,1)));
yMin = min(min(x1(:,2)),min(x2(:,2)));
zMin = min(min(x1(:,3)),min(x2(:,3)));

%Output

fig1 = figure(1);
set(fig1, 'position', [200,200,1000,700]);
plot(t1,x1(:,1),'-b',t2,x2(:,1),'--k','LineWidth',1.5);
xlabel('$t$','FontSize',16,'Interpreter','LaTex');
ylabel('$u$-value','FontSize',16,'Interpreter','LaTex');
% legend({'Initial condition $(1,1,1)^T$',...
% 	'Initial condition $(0.9999998,1.0000002,1.0000002)^T$'},...
% 	'FontSize',13,'Interpreter','LaTex');
axis([0 50 -30 30]);

print('-depsc2', 'LorenzU01.eps','-b0'); 
print('-dpdf', 'LorenzU01.pdf','-b0');


fig2 = figure(2);
set(fig2, 'position', [200,200,1200,700]);
sub1 = subplot(1,2,1);
plot3(x1(:,1),x1(:,2),x1(:,3),'-b','LineWidth',1.2);
xlabel('$u$','FontSize',16,'Interpreter','LaTex');
ylabel('$y$','FontSize',16,'Interpreter','LaTex');
zlabel('$z$','FontSize',16,'Interpreter','LaTex');
%title('Initial condition  $(1,1,1)^T$','FontSize',14,...
%        'Interpreter','LaTex');
ax = gca;
axis equal;
grid on;
axis([xMin, xMax, yMin, yMax, zMin, zMax]);
set(ax,'Ydir','reverse');

sub2 = subplot(1,2,2);
plot3(x2(:,1),x2(:,2),x2(:,3),'-k','LineWidth',1.2);
xlabel('$u$','FontSize',16,'Interpreter','LaTex');
ylabel('$y$','FontSize',16,'Interpreter','LaTex');
zlabel('$z$','FontSize',16,'Interpreter','LaTex');
% title('Initial condition $(0.9999998,1.0000002,1.0000002)^T$',...
%		'FontSize',14,'Interpreter','LaTex');
ax = gca;
axis equal;
grid on;
axis([xMin, xMax, yMin, yMax, zMin, zMax]);
set(ax,'Ydir','reverse');

print('-depsc2', 'LorenzTrajectory01.eps','-b0'); 
print('-dpdf', 'LorenzTrajectory01.pdf','-b0');

