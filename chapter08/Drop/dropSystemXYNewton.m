%% IVP for drop
%% theta = t 
%% System for x and y

clear all; close all; 

ltype = {'b-','r--','m-.','k:'};             % for Plot

gamma = 1.;
rho= 1.;
g  = 10.;
p0 = 2.;
a = sqrt(rho*g/gamma); 
p = p0/(a*gamma);

% hydrophobic
thetac=2*pi/3;

% ﻿hydrophilic
% thetac=pi/4;

KN=100;
N = 1000; 
h = thetac /(N-1) ; 

x = ones(N,1) ;
y = ones(N,1) ;

x(1)=0.0; 
y(1)=0.0; 

% Theta loop
for j=2:N
%theta    
    tj = (j-1)*h;
    
% Newton solve  
    xk = x(j-1);
    yk = y(j-1);
    
    if j == 2
        xk = 10.0*h; 
        yk = xk * tan(h);
    end
for k = 1:KN

    Ek = xk - h * dropfXY(xk,yk,tj,p) * cos(tj) - x(j-1);
    dEkdx = 1 - h * dropdfdxXY(xk,yk,tj,p) * cos(tj);
    dEkdy = - h * dropdfdyXY(xk,yk,tj,p) * cos(tj);
    
    Fk = yk - h * dropfXY(xk,yk,tj,p) * sin(tj) - y(j-1); 
    dFkdx =  - h * dropdfdxXY(xk,yk,tj,p) * sin(tj);
    dFkdy = 1 - h * dropdfdyXY(xk,yk,tj,p) * sin(tj);
    
    % Jacobian 
    JEF = [ dEkdx , dEkdy ; dFkdx , dFkdy ]  ; 
    
    sk = [ xk ; yk ]; 
    
    sk = sk - inv(JEF) * [Ek ; Fk] ; 
    
    xk = sk(1) ;
    yk = sk(2) ;   
    
    if (Ek^2 + Fk^2) < 10e-15 
        
        break
    end
end    
    x(j) = xk  ;
    y(j) = yk  ;

end


ybottom = - y(N);

% plot the solution in 'a' units
figure(1)

plot(x,-y,ltype{1},'Linewidth',1); hold on; 
plot(-x,-y,ltype{1},'Linewidth',1); hold on; 


yline(ybottom,ltype{2},'Linewidth',1);

print('-depsc2', 'drop.eps','-b0'); 
 print('-dpdf', 'drop.pdf','-b0');

