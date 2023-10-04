
close all; 
clear all;

N = 30000;
% h = 0.0001 ;
K = 1;  
m0= 1; 
c = 1; 
%angular frequency
omega0 = sqrt(K/m0);
T0 = 2 * pi / omega0; 

h = 2 * T0 / N ;
%initial conditions
Xim1 = 0; 
Uim1 = 0.98; 

for i = 1:N 


T(i) = i*h ;
X(:,i) = Xim1 + h * Uim1;

Ai = (h / m0) * (1 - Uim1^2 / c^2 ) * (- K * Xim1) ;

U(:,i) = (c^2 * Uim1 + c * Ai * sqrt(c^2 - Uim1^2 + Ai^2)) / (c^2 + Ai^2) ;

Xim1 = X(:,i); 
Uim1 = U(:,i); 

end; 


%  
ltype = {'b-','r--','m-.'}; 

  figure(1)
  
  grid on; hold on; 
   plot (T,X, ltype{1},'Linewidth',2); 
   plot (T,U, ltype{2},'Linewidth',2);
   

   print('-depsc2', 'oscillatorSR01.eps','-b0'); 
   print('-dpdf', 'oscillatorSR01.pdf','-b0');


   
%  xlabel('Time t'); ylabel('H(q,p) -  H(q_0,p_0)');
% plot (T,error_L) ; grid on;
%  xlabel('Time t');  ylabel('|| L(q,p) -  L(q_0,p_0)||');
% plot (T,error_RL) ;grid on;
%   xlabel('Time t'); ylabel('|| A(q,p) -  A(q_0,p_0)||');

