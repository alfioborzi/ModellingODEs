% Bargaining - Nash versus Pareto
% Our matrix Riccati differential equation:
% - dX/dt = A'X + XA + X B R^-1 B' X - Q ; X(T)=XT
% We reverse time : tau = T - t 
% Our transformed Riccati equation
% dX1/d tau = A' X1 + X1 A + X1 E1 X1 + X1 E2 X2 + X2 E2 X1 - Q1; X1(0)=XT1
% dX2/d tau = A' X2 + X2 A + X2 E2 X2 + X1 E1 X2 + X2 E1 X1 - Q2; X2(0)=XT2
%
% We consider two players 

close all; clear all;

a1 = inline('lambda','lambda'); 
a2 = inline('1-lambda','lambda');
 
A = 1; 

B1 = 1; 
Q1 = 2;
R1 = 1;
XT1 = -1;
D1 = - XT1; 
 
B2 = 2; 
Q2 = 1;
R2 = 2;
XT2 = -2;
D2 = - XT2;

% Pareto
B = [B1 B2];
% Assemble XT for Nash
XT = [XT1; XT2]; 


% Notice that R1 and R2 are elements of matrix
R1P = [R1 0; 0 0];
R2P = [0 0; 0 R2];

% Final time, t0=0 initial time
T  = 1.0; t0=0; 
% Let S0 be the initial state 
S0 = 2;


% consider lambda = m*Dlambda, m=1,...,Nm, Dlambda=1/Nm
Nm = 30;
Nm1= Nm-1;
Dlambda = 1/Nm;
lambda0 = 0;

% Nash product
 NBP  = 0;
 NBPV = zeros(Nm1);

fileID = fopen('lambdaPO.txt','w');

% Time grid from t0 to T
dt = 0.0001;
Time = t0:dt:T;
iTime=size(Time,2);

%----------- determine the Nash point --------

% Nash: For efficiency precompute E = B R^-1 B'
R1i = 1/R1;
E1 = B1*R1i*B1;

R2i = 1/R2;
E2 = B2*R2i*B2;

% ode4 for Nash
[X] = ode5(@(t,X)mRiccatiNash(t, X, A, E1, E2, Q1, Q2), Time, XT);

% Assembling the value functions V(S0,t0)

%---------------------------------------
% Nash
% assemble A + E * X 

nX  =size(X,2);
nXm =nX/2; 

ACN=zeros(iTime,1); 
for i=1:1:iTime
X1 = reshape(X(iTime-i+1,1:nXm),size(A));
X2 = reshape(X(iTime-i+1,nXm+1:nX),size(A));
AMN = A + E1*X1 + E2*X2;
ACN(i,:) = AMN(:); 
end

% State with initial condition y0 at t0 and Nash controls
y0=S0;
[y] = ode5(@(t,y)statedyn(t, y, ACN, A, t0, dt), Time, y0);
%
% figure(6)
% plot(y)
% title('State y - Nash solution')

% Nash controls corresponding to y 
un = zeros(iTime,2);
for i=1:1:iTime
% XX1 = X(iTime-i+1,1);
% XX2 = X(iTime-i+1,2);
XX1 = reshape(X(iTime-i+1,1:nXm),size(A));
XX2 = reshape(X(iTime-i+1,nXm+1:nX),size(A));
YS  = y(i);
u1(i,:) = R1i*B1.'*XX1*YS;
u2(i,:) = R2i*B2.'*XX2*YS;
un(i,:)=[u1(i,:) u2(i,:)];
end
%
% figure(7)
% plot(un)
% title('Control u - Nash solution')

% Nash: compute the value functions - cost functionals
% for initial state S0 at t0 
J1N = 0; 
for i=1:1:iTime-1
J1N = J1N + (y(i)*Q1*y(i) + un(i,1)*R1*un(i,1));
end
J1N = dt * J1N + y(iTime)*D1*y(iTime);


J2N = 0; 
for i=1:1:iTime-1
J2N = J2N + (y(i)*Q2*y(i) + un(i,2)*R2*un(i,2));
end
J2N = dt * J2N + y(iTime)*D2*y(iTime);

% (J1N,J2N) is the Nash point

%---------------------------------------



% ---- determine the Pareto front ----

for m=2:Nm
lambda = (m-1)*Dlambda + lambda0;

% Pareto averaging
Q = a1(lambda)*Q1+a2(lambda)*Q2;
R = a1(lambda)*R1P+a2(lambda)*R2P;
YT = a1(lambda)*XT1+a2(lambda)*XT2;

% Pareto: For efficiency precompute E = B R^-1 B'
Ri = inv(R);
E = B*Ri*B.';

% ode4 for Pareto; uses fix time step size 
[Y] = ode5(@(t,Y)mRiccati(t, Y, A, E, Q), Time, YT);
%
%---------------------------------------
% Pareto
% assemble Delta = A + E * Y 
%

ACP=zeros(size(Y));
for i=1:1:iTime
YY = reshape(Y(iTime-i+1,:), size(A));
AMP = A + E*YY;
ACP(i,:) = AMP(:); 
end
% State with initial condition y0 at t0 and Pareto control
y0=S0;
[y] = ode5(@(t,y)statedyn(t, y, ACP, A, t0, dt), Time, y0);
%

% Pareto controls corresponding to y 
up = zeros(iTime,2);
for i=1:1:iTime
YY = Y(iTime-i+1);
YS = y(i);
up(i,:) = Ri*B.'*YY*YS;
end
%


% Pareto: compute the cost functionals
% for initial state S0 at t0 
J1P = 0; 
for i=1:1:iTime-1
J1P = J1P + (y(i)*Q1*y(i) + up(i,1)*R1*up(i,1));
end
J1P = dt * J1P + y(iTime)*D1*y(iTime);

J2P = 0; 
for i=1:1:iTime-1
J2P = J2P + (y(i)*Q2*y(i) + up(i,2)*R2*up(i,2));
end
J2P = dt * J2P + y(iTime)*D2*y(iTime);

% (J1P,J2P) are the points on the Pareto front

% check bargaining
 
NBP=(J1P-J1N)*(J2P-J2N);

fprintf('lambda=%8.4f Nash distance=%8.4f \n',lambda,NBP);


NBPV(m)=NBP;
 
fprintf(fileID,'%8.4f %8.4f %8.4f %8.4f %8.4f  %8.4f \n',lambda,J1P,J2P,J1N,J2N,NBP);


end  

% The Nash point
NE = [J1N,J2N]

% The PO point maximizing NBP
[M,I] = max(NBPV);
vPO=M(1);
mPO=I(1);
lambdaNB = (mPO-1)*Dlambda + lambda0 
% Pareto averaging
Q = a1(lambda)*Q1+a2(lambda)*Q2;
R = a1(lambda)*R1P+a2(lambda)*R2P;
YT = a1(lambda)*XT1+a2(lambda)*XT2;

% Pareto: For efficiency precompute E = B R^-1 B'
Ri = inv(R);
E = B*Ri*B.';

% ode4 for Pareto; uses fix time step size 
[Y] = ode5(@(t,Y)mRiccati(t, Y, A, E, Q), Time, YT);
ACP=zeros(size(Y));
for i=1:1:iTime
YY = reshape(Y(iTime-i+1,:), size(A));
AMP = A + E*YY;
ACP(i,:) = AMP(:); 
end
% State with initial condition y0 at t0 and Pareto control
y0=S0;
[y] = ode5(@(t,y)statedyn(t, y, ACP, A, t0, dt), Time, y0);
%
% figure(4)
% plot(y)
% title('State y - Pareto solution at baragining point')

% Pareto controls corresponding to y 
up = zeros(iTime,2);
for i=1:1:iTime
YY = Y(iTime-i+1);
YS = y(i);
up(i,:) = Ri*B.'*YY*YS;
end
% figure(5)
% plot(up)
% title('Control u - Pareto solution at bargaining point')

% Pareto at bargaining point: compute the cost functionals
%
J1P = 0; 
for i=1:1:iTime-1
J1P = J1P + (y(i)*Q1*y(i) + up(i,1)*R1*up(i,1));
end
J1P = dt * J1P + y(iTime)*D1*y(iTime);

J2P = 0; 
for i=1:1:iTime-1
J2P = J2P + (y(i)*Q2*y(i) + up(i,2)*R2*up(i,2));
end
J2P = dt * J2P + y(iTime)*D2*y(iTime);

% (J1P,J2P) the PO point at bargaining
PB=[J1P,J2P]
%---------

 fclose(fileID);
 
 
 figure(1)
 load lambdaPO.txt
 [n,p] = size(lambdaPO);
 plot(lambdaPO(:,2),lambdaPO(:,3),'LineWidth',3,'Color','b')
 hold on
 plot(J1N,J2N,'k*');
 plot(J1P,J2P,'ob','MarkerSize',10);

print('-depsc2', 'POfrontODE01.eps','-b0'); 

print('-dpdf', 'POfrontODE01.pdf','-b0');

 


%--------- 

% End Main 

%--------- Set up ODEs -----------
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

function dXdt = mRiccati(t, X, A, E, Q)
% Our transformed Riccati equation
% dX/d tau = A'X + XA + X B R^-1 B' X - Q ; 
dXdt = A.'*X + X*A + X*E*X - Q; %Determine derivative
dXdt = dXdt(:); %Convert from "n"-by-"n" to "n^2"-by-1
end

%--------- 

function dXdt = mRiccatiNash(t, X, A, E1, E2, Q1, Q2)
% Our transformed coupled Riccati equations
% dX1/d tau = A' X1 + X1 A + X1 E1 X1 + X1 E2 X2 + X2 E2 X1 - Q1 ; 
% dX2/d tau = A' X2 + X2 A + X2 E2 X2 + X1 E1 X2 + X2 E1 X1 - Q2 ; 
X1 = X(1);
X2 = X(2);
dX1dt = A.'*X1 + X1*A + X1*E1*X1 + X1*E2*X2 + X2*E2*X1 - Q1; %Determine derivative
dX2dt = A.'*X2 + X2*A + X2*E2*X2 + X1*E1*X2 + X2*E1*X1 - Q2; %Determine derivative
dXdt = [dX1dt(:); dX2dt(:)] ; %Convert from 2 "n"-by-"n" to "2 n^2"-by-1
end

%--------- 

function dZdt = mRiccatiPareto(t, Z, AC, AF, KC, E, Q, t0, dt)
% auxiliary Riccati equation: Ehtamo (41a)
% dX/d tau = A'X + XA + X B R^-1 B' X - Q ; 
% time step - back
jf = ((t-t0)/dt);
j = uint32(jf)+1;
% A(t)
AT = AC(j,:);
D = reshape(AT, size(AF)); %Convert from "n^2"-by-1 to "n"-by-"n"
KT = KC(j,:);
K = reshape(KT, size(AF)); %Convert from "n^2"-by-1 to "n"-by-"n"

dZdt = D.'*Z + Z*D - K*E*K - Q; %Determine derivative
dZdt = dZdt(:); %Convert from "n"-by-"n" to "n^2"-by-1
end
