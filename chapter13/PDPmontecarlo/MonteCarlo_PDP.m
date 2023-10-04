%function MonteCarlo_PDP(soluz_processo,paramsA,x0,MarkovP,T,Ntraiet,seq,bname)
%
% MonteCarlo_PDP(soluz_processo,paramsA,x0,MarkovP,T)
% visualize trajectories of a PDP process 
%
% soluz_processo: file with analytical sol 
% paramsA: parameters [gam1 gam2 ; W1 W2]
% x0: initial state
% MarkovP: matrix of transitions 
% seq: save 1 in seq pictures
% bname: name of picture files, if empty not saved

% t(1)=0;
% x(1)=x0;

close all;
clear all; 

ltype = {'b-','r--','m-.'}; 

x0=0;    % inital state
T=10;    % final time
Ntraiet = 3;  % number of realisations 


% parameter of deterministic dynamics 
%   -gam*x+W  and -gam*x-W
param.gam = [ 4 4 ];  
param.W   = [ -2 2];

% comp domain [a,b]
a=param.W(1)/param.gam(1);
b=param.W(2)/param.gam(2);

% parameter Markovian process
q=[0 1 ; 1 0]; 
lamb=[2 2]';

stati=length(lamb);   % number of states


% cumulative of transition matrix 
 cumul = zeros(stati);
   for s=1:stati
	 cumul(1,s)=q(1,s);
     for j=2:stati
         cumul(j,s)=cumul(j-1,s)+q(j,s);
     end
   end

n=0;  
   

while(Ntraiet)

t1=0;
x1=0;

s=floor(rand*stati)+1;


%hold off;
while(t1<T)
    dt=genera_Poisson(lamb(s));
   
    [x2,t2]=traiet_filter(x1,t1,s,dt,param);
        
    % sort a new dyn state
    r=rand;
	  
    % find the corresponding matrix element
    for j=1:stati
         if r < cumul(j,s)
            s=j;
	    break;
         end
    end
      
        plot(t2,x2,ltype{Ntraiet},'Linewidth',2);
     
        axis([0 T a b]);
            
        pause(0.04);    
     
     n=n+1;
     
     t1=t2(end);
     x1=x2(end);    
     hold on;
%     pause(0.5);
end

Ntraiet = Ntraiet-1;

end 

ybottom = 0;
yline(ybottom,'b','Linewidth',1);


print('-depsc2', 'PDP01.eps','-b0'); 
print('-dpdf', 'PDP01.pdf','-b0');


hold off;
