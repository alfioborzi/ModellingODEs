
close all;
clear all;

ltype = {'b-','r--','m-.','k:'};             % for Plot


f = @(t,y) [y(2); -sin(y(1)) - y(2)];


y1 = linspace(-7,7,20);
y2 = linspace(-2,2,20);

% creates two matrices one for all the x-values on the grid, and one for
% all the y-values on the grid. Note that x and y are matrices of the same
% size and shape, in this case 20 rows and 20 columns
[x,y] = meshgrid(y1,y2);


u = zeros(size(x));
v = zeros(size(x));

% we can use a single loop over each element to compute the derivatives at
% each point (y1, y2)
t=0; % we want the derivatives at each point at t=0, i.e. the starting time
for i = 1:numel(x)
    Yprime = f(t,[x(i); y(i)]);
    u(i) = Yprime(1);
    v(i) = Yprime(2);
end

% This way the arrows are normalised 
for i = 1:numel(x)
Vmod = sqrt(u(i)^2 + v(i)^2);
u(i) = 0.1*u(i)/Vmod;
v(i) = 0.1*v(i)/Vmod;
end

figure(1)
quiver(x,y,u,v,'r'); 
xlabel('y_1')
ylabel('y_2')
axis tight equal;

%figure(2)
 hold on
for y1i = -7:1:7
for y2i = -2:1:2
    [ts,ys] = ode45(f,[0,5],[y1i;y2i]);
    plot(ys(:,1),ys(:,2), ltype{1},'Linewidth',1);
    hold on;
end
end
xL = xlim;
yL = ylim;
line([0 0], yL);  %y-axis
line(xL, [0 0]);  %x-axis

hold off

print('-depsc2', 'trajNL01.eps','-b0'); 
print('-dpdf', 'trajNL01.pdf','-b0');

