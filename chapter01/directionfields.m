
close all;
clear all;

% dx/dt = x^2-3xy+y
% dy/dt = -5x+sin(yx)
% 

[x,y] = meshgrid(-2:0.2:2);
dx = x.^2-3*x.*y+y;
dy = -5*x+sin(x.*y);
r = ( dx.^2 + dy.^2 ).^0.5;
px = dx./r;
py = dy./r;


figH =figure; 
[x,y] = meshgrid(-4:.5:4, -4:0.5:4);
dy = x - sin(y); 
dx = ones(size(dy)); 
r = ( dx.^2 + dy.^2 ).^0.5;
px = dx./r;
py = dy./r;
figure(1)
quiver(x,y,px,py);
%print('-depsc2', 'exDirField01.eps'); 

print('-dpdf', 'exDirField01.pdf','-r0');


