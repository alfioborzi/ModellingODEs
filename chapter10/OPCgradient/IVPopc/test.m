

clear all 
close all

% Define Parameter of the Control Problem
alpha=2;
beta=2;
nu=10^-2;
y0=1;
yT=2;
a=-3;
B=1;
t0=0;
T=1.5;

% Choose number of discrete time steps
N=100;

h=T/N;


% Choose tolerance and maximal number of iterations
epsilon=10^-2;
maxiter=50;


% Solves the optimal control problems
[y1,p1,u1,normdJ1,J1,counter1] = optimalexpliciteuler (a,B,t0,T,N,alpha,beta,nu,y0,yT,epsilon,maxiter);


% plot
subplot(1,2,1)
hold on
x1=[0:h:T];
%controlled trajectory
plot(x1,y1,'r')
%desired trajectory
plot(x1,yd(x1),'g')
% initial value
plot(0,y0,'*')
% desired terminal value
plot(T,yT,'*')
title ('y');


% subplot(2,3,2)
% hold on
% plot(x1,p1,'r')
% 
% title ('p');
% 
% subplot(2,3,3)
% hold on
% plot(x1,u1,'r')
% 
% title ('u');

subplot(1,2,2)
hold on
x2=[0:1:counter1-1];
plot(x2,normdJ1,'r')
title ('Norm of dJ');

% subplot(2,3,5)
% hold on
% x2=[0:1:counter1];
% plot(x2,J1,'r')
% title ('J');

hold off


