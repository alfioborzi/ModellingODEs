
close all; 
clear all;

N = 3500;
h = 0.01 ;
K = 1;  
%initial conditions
q0 = [1.0,0.5,0.0]'; 
p0 = [0.0,1.0,0.5]'; 
r0 = sqrt( q0'*q0) ;
q = q0 ;
p = p0 ;
r = sqrt( q'*q) ;
%first energy. If H0 <0 then ellipses
H0= ( (p'*p/2)-(K/r)) ;
%for angular momentum 
m1 = q(2)*p(3) - q(3)*p(2);
m2 = q(3)*p(1) - q(1)*p(3); 
m3 = q(1)*p(2) - q(2)*p(1);
%first angular momentum
L0 =[m1,m2,m3]';
%first Runge-Lenz vector
RL0 = q *(H0 + p'*p/2 ) - p *(q'*p) ; 
% computational unit

for i = 1:N 
q1 = q ; 
p1 = p; 
q2 = q1 ;
p2 = p1;
% implicit solve
for k = 1:100 
p_mid = (p+p1)/2;
q2= q+h*p_mid;
q_mid = (q+q1)/2;
r = sqrt( q_mid'*q_mid);
p2 = p - h*K*q_mid/(r^3) ; 
q1= q2;
p1 = p2;
end; 
Q(:,i) = q1;
P(:,i) = p1;
T(i) = i*h ;
q = q1; p = p1; 
H(i) = ( (p'*p/2)-(K/sqrt(q'*q))); 
error_H(i) = (H(i) - H0) ; 

%for angular momentum

m1 = q(2)*p(3) - q(3)*p(2);
m2 = q(3)*p(1) - q(1)*p(3);
m3 = q(1)*p(2) - q(2)*p(1);
L(:,i)= [m1,m2,m3]' ;
error_L(:,i) = sqrt( (L(:,i) - L0)'*(L(:,i)-L0)) ; 
RL(:,i)=q*(H(i)+p'*p/2)-p*(q'*p) ; 
error_RL(:,i) = sqrt( (RL(:,i)- RL0)'*(RL(:,i)- RL0)) ; 
end; 


figure(1)
plot3 (Q(1,:),Q(2,:),Q(3,:) ) ; grid on; hold on; 
plot3 (0,0,0, 'r*' ) ; 
  xlabel('q_1');
  ylabel('q_2');
  zlabel('q_3');
  
print('-depsc2', 'Kepler01.eps','-b0'); 
print('-dpdf', 'Kepler01.pdf','-b0');

%  
ltype = {'b-','r--','m-.'}; 

  figure(2)
  
  grid on; hold on; 
   plot (T,error_H, ltype{1},'Linewidth',2); 
   plot (T,error_L, ltype{2},'Linewidth',2);
   plot (T,error_RL, ltype{3},'Linewidth',2);

   print('-depsc2', 'KeplerErrMP01.eps','-b0'); 
   print('-dpdf', 'KeplerErrMP01.pdf','-b0');


   

