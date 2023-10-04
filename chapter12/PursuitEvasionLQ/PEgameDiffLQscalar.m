% 
% 
% Pursuit-Evasion in a LQ diff game setting
%

close all; clear all;

% for reference to the book
% 1 -> e , 2 -> p

ltype = {'b-','r--','m-.','k:'};             % for Plot


% parameters

c_e = 1.0e+0 ; c_ep = 1.0e+6; sg_e = 4.0e+0; sg_e2 = sg_e * sg_e; 
c_p = 4.0e+0 ; c_pe = 1.0e+6; sg_p = 1.0e+0; sg_p2 = sg_p * sg_p; 

% Final time, t0=0 initial time
T  = 1.0; t0=0; 

% Define position and velocities in 2D : [r1,r2,v1,v2]
% Initial position and velocity of the pursuer
yp = [1 ; 1 ; 0.5 ; 0.5] ;

% Initial position and velocity of the evader
ye = [1.5 ; 0.75 ; 0.125 ; 1.5];

% 2 x 2 identity
I2 = eye(2); 
% 2 x 2 zero
I0 = zeros(2,2);

% dynamics
A = [I0 I2; I0 I0];

[n,m] = size(A);

% y0 initial state of relative position and velocity: y0=ye-yp
y0 = ye - yp ;


% e functional and control
B1 = [I0; I2];
M1 = [I0 I0; I0 I0];
N11 = (1/c_e) * I2;
N12 = (1/c_ep) * I2;
D1 = - sg_e2 * [I2 I0; I0 I0];

% p functional and control
 
B2 = [I0; -I2];
M2 = [I0 I0; I0 I0];
N21 = (1/c_pe) * I2;
N22 = (1/c_p) * I2;
D2 = sg_p2 * [I2 I0; I0 I0];

%----------- determine the Nash point --------

% Time grid from t0 to T
N  = 1000;
iTime=N+1;
dt = T/N;
Time = t0:dt:T;

% Terminal condition for the scalar Riccati problem
wT=[c_e * sg_e2 ; - c_p * sg_p2 ]; 

c1=c_e/c_ep; 
c2=c_p/c_pe;

Nw=1000; 

Tw = T^3/3. ;
dtw = Tw/(Nw*N) ;
Timew = t0:dtw:Tw;

[w] = ode5(@(t,w)peRiccatiScalar(t,w,c1,c2), Timew, wT);


for i=1:1:iTime
    irev = Nw*(iTime-i)+1;
w1(i) = w(irev,1)/c_e;
w2(i) = w(irev,2)/c_p;
end

%---------------------------------------
% Controlled system
% assemble A + B * N^-1 B' * Q 
% define C = B * N^-1 B'

C1 = B1*inv(N11)*B1';
C2 = B2*inv(N22)*B2';

% compute trajectories
for i=1:1:iTime

tt = T - Time(i);

Y1 = w1(i) * [ I2  tt*I2 ;  tt*I2  tt^2 * I2 ];
Y2 = w2(i) * [ I2  tt*I2 ;  tt*I2  tt^2 * I2 ];

% relative motion
AMN = A + C1*Y1 + C2*Y2;
ACN(i,:) = AMN(:); 
end

% State of relative position and velocity 
% with initial condition y0 at t0 and Nash controls

[y] = ode5(@(t,y)statedyn(t, y, ACN, A, t0, dt), Time, y0);


 figure(1)
 plot([t0:dt:T],y(:,1),ltype{1},'Linewidth',1); hold on; 
 plot([t0:dt:T],y(:,2),ltype{2},'Linewidth',1); hold on; 
 plot([t0:dt:T],y(:,3),ltype{3},'Linewidth',1); hold on; 
 plot([t0:dt:T],y(:,4),ltype{4},'Linewidth',1); hold off; 
 
 legend({'$r_1$','$r_2$',...
	'$v_1$','$v_2$'},'Interpreter','Latex','FontSize',12)
 
% title('Relative state y = (r,v) - Nash Eq. solution')
 
print('-depsc2', 'reLQstate.eps','-b0'); 
% print('-dpdf', 'reLQstate.pdf','-b0');



% Nash controls corresponding to y 

for i=1:1:iTime
    
%tt = Time(iTime-i+1);

tt = T - Time(i);


%w1 = w(iTime-i+1,1)/c_e;
%w2 = w(iTime-i+1,2)/c_p;

Y1 = w1(i) * [ I2  tt*I2 ;  tt*I2  tt^2 * I2 ];
Y2 = w2(i) * [ I2  tt*I2 ;  tt*I2  tt^2 * I2 ];

YS  = y(i,:)';

u1(i,:) = inv(N11)*B1.'*Y1*YS;
u2(i,:) = inv(N22)*B2.'*Y2*YS;

end
%
 figure(2)
 
 plot([t0:dt:T],u1(:,1),ltype{1},'Linewidth',1); hold on; 
 plot([t0:dt:T],u1(:,2),ltype{2},'Linewidth',1); hold on; 
 plot([t0:dt:T],u2(:,1),ltype{3},'Linewidth',1); hold on; 
 plot([t0:dt:T],u2(:,2),ltype{4},'Linewidth',1); hold off; 
 
 legend({'$a_e (1)$','$a_e (2)$',...
	'$a_p (1)$','$a_p (2)$'},'Interpreter','Latex','FontSize',12)

 
% title('Control  - Nash Eq. solution')
 
print('-depsc2', 'peLQcontrol.eps','-b0'); 
% print('-dpdf', 'peLQcontrol.pdf','-b0');


% State of position and velocity of the evader 
% with initial condition ye at t0 and acceleration = Nash control

 [ze] = ode5(@(t,ze)epdyn(t, ze, A, B1, u1, t0, dt), Time, ye);
%
% State of position and velocity of the pursuer 
% with initial condition yp at t0 and acceleration = Nash control

 [zp] = ode5(@(t,zp)epdyn(t, zp, A, B1, u2, t0, dt), Time, yp);
%


figure(4) 

 plot(ze(:,1),ze(:,2),ltype{1},'Linewidth',1); hold on; 
 plot(zp(:,1),zp(:,2),ltype{2},'Linewidth',1); hold off; 
 
 
 legend({'$r_e$','$r_p$'},'Interpreter','Latex','FontSize',12,'Location','southeast')
 
 
print('-depsc2', 'peLQstate.eps','-b0'); 
% print('-dpdf', 'peLQstate.pdf','-b0');


figure(5)

plot([t0:dtw:Tw],w(:,1),ltype{1},'Linewidth',1); hold on; 
plot([t0:dtw:Tw],w(:,2),ltype{2},'Linewidth',1); hold off; 

 legend({'$\eta$','$\pi$'},'Interpreter','Latex','FontSize',12,'Location','southeast')

 print('-depsc2', 'peLQetapi.eps','-b0'); 
% print('-dpdf', 'peLQetapi.pdf','-b0');

%--------- --------- --------- -------- --------- 

% End Main 

%--------- Set up ODE system for e/p motion -----------

function dydt = epdyn(t, y, AF, BF, u, t0, dt)
% State equation
% dy/dt = A(t) * y + B * u; 
% time step

jf = ((t-t0)/dt);
j = uint32(jf)+1;
% A(t)
ut = u(j,:);
dydt = AF*y + BF*ut' ;                %Determine derivative
end

%--------- Set up ODE system for relative motion -----------
function dydt = statedyn(t, y, AC, AF, t0, dt)
% State equation
% dy/dt = A(t) * y ; 
% time step
jf = ((t-t0)/dt);
j = uint32(jf)+1;
% A(t)
AT = AC(j,:);
AD = reshape(AT, size(AF)); %Convert from "n^2"-by-1 to "n"-by-"n"
dydt = AD*y;                %Determine derivative
end

%--------- 

function dwdt = peRiccatiScalar(t,w,c1,c2)
%c1=c_e/c_ep; 
%c2=c_p/c_pe;
dw1dt = w(1)^2 + 2*w(1)*w(2) - c1 * w(2)^2 ; 
dw2dt = w(2)^2 + 2*w(1)*w(2) - c2 * w(1)^2 ; 
dwdt = [dw1dt ; dw2dt] ; 
end



