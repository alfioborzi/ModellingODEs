function [ ] = SQHquantum( )
% SQH Scheme
close all; 
clear all;

ltype = {'b-','r--','m-.','b-*','r:','m-x'};

OCP=struct('dt',10^-5,'T',0.01,'umin',-60,'umax',60,'nu',10^-6,'beta',10^-2);
y0=[0;0;1;0;0;1];
y0=y0/norm(y0);
A=[0,-1,0,0,0,0;1,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,1,0;0,0,0,-1,0,0;0,0,0,0,0,0]*2*pi*483;
B=[0,0,0,0,0,0;0,0,-1,0,0,0;0,1,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,-1;0,0,0,0,1,0]*2*pi;

%Tolerance for convergence
kappa=10^-15;                
%Algorithm paramters
eta=10^-9;
zeta=0.8;
sigma=2.0;
kmax=1000;        %Maximum iteration number
epsilon=10.;      %Initial guess for epsilon

yd=[0;1;0;0;1;0];
yd=yd/norm(yd);
Nt=round(OCP.T/OCP.dt);

u=0*ones(1,Nt+1);
u_old=u;

y=forward(A,B,y0,u,OCP);
y_old=y;

p=backward(A,B,y(:,end),yd,u,OCP);

J_plot(1,1)=get_J(y(:,end),u,yd,OCP);
J_plot(1,1);
Jk(1)=J_plot(1,1);

count_updates=1;
tic
for k=1:kmax
    for i=1:Nt+1
        %u(1,i)=argminH(y(:,i),u(i),p(:,i),A,B,epsi,OCP );

        uz(1)=min(max(0,(2*epsilon*u_old(1,i)-(p(:,i)')*B*y(:,i)-OCP.beta)/(2*epsilon+OCP.nu)),OCP.umax);
        uz(2)=min(max(OCP.umin,(2*epsilon*u_old(1,i)-(p(:,i)')*B*y(:,i)+OCP.beta)/(2*epsilon+OCP.nu)),0);
        H(1)=(OCP.nu/2)*uz(1)^2 + OCP.beta*abs(uz(1)) + (p(:,i)')*(A+uz(1)*B)*y(:,i) + epsilon*(uz(1)-u_old(1,i))^2;
        H(2)=(OCP.nu/2)*uz(2)^2 + OCP.beta*abs(uz(2)) + (p(:,i)')*(A+uz(2)*B)*y(:,i) + epsilon*(uz(2)-u_old(1,i))^2;
        [~,pos]=min(H);
        u(1,i)=uz(pos);

    end 
  
   du=sum((u-u_old).^2)*OCP.dt;
   
   y=forward(A,B,y0,u,OCP);
   J_int=get_J(y(:,end),u,yd,OCP);
   if(J_int-Jk(count_updates)>-eta*du)
       u=u_old;
       y=y_old;
       epsilon=epsilon*sigma;                            
   else
       count_updates=count_updates+1;
       p=backward(A,B,y(:,end),yd,u,OCP);
       epsilon=epsilon*zeta;                        
       Jk(count_updates)=J_int;
       u_old=u;
       y_old=y;
   end
   if(du<kappa)
        u=u_old;
        fprintf('converged\n')
        break;
   end 
end 
toc
figure(1)
plot(0:OCP.dt/(OCP.T):1,u/100,ltype{1},'Linewidth',2)
axis([-inf inf -1 1])
xlabel('x/T')
ylabel('u/100')
print('-depsc2', 'quantumOPC02.eps','-b0'); 
print('-dpdf', 'quantumOPC02.pdf','-b0');

figure(2)
plot(0:OCP.dt/(OCP.T):1,y(1,:),ltype{1},'Linewidth',2), hold on; 
plot(0:OCP.dt/(OCP.T):1,y(2,:),ltype{2},'Linewidth',2),
plot(0:OCP.dt/(OCP.T):1,y(3,:),ltype{3},'Linewidth',2),
%plot(0:OCP.dt/(OCP.T):1,y(4,:),ltype{4},'Linewidth',1),
%plot(0:OCP.dt/(OCP.T):1,y(5,:),ltype{5},'Linewidth',1),
%plot(0:OCP.dt/(OCP.T):1,y(6,:),ltype{6},'Linewidth',1),
axis([-inf inf -1 1])
xlabel('x/T')
ylabel('y')
legend({'$y_1$','$y_2$','$y_3$','$y_4$','$y_5$',...
	'$y_6$'},'Interpreter','Latex','FontSize',12)

print('-depsc2', 'trajectory13OPC02.eps','-b0'); 
print('-dpdf', 'trajectory13OPC02.pdf','-b0');

figure(3)
plot(0:OCP.dt/(OCP.T):1,y(4,:),ltype{1},'Linewidth',2), hold on; 
plot(0:OCP.dt/(OCP.T):1,y(5,:),ltype{2},'Linewidth',2),
plot(0:OCP.dt/(OCP.T):1,y(6,:),ltype{3},'Linewidth',2),
%plot(0:OCP.dt/(OCP.T):1,y(4,:),ltype{4},'Linewidth',1),
%plot(0:OCP.dt/(OCP.T):1,y(5,:),ltype{5},'Linewidth',1),
%plot(0:OCP.dt/(OCP.T):1,y(6,:),ltype{6},'Linewidth',1),
axis([-inf inf -1 1]);
xlabel('x/T')
ylabel('y')
legend({'$y_4$','$y_5$',...
	'$y_6$'},'Interpreter','Latex','FontSize',12)

print('-depsc2', 'trajectory13OPC02.eps','-b0'); 
print('-dpdf', 'trajectory13OPC02.pdf','-b0');


figure(4)
plot(Jk,ltype{2},'Linewidth',2)
xlabel('SQH iterations')
ylabel('J')

C=0.5*((y(:,end)-yd)'*(y(:,end)-yd));

numeric_optimality(u,y0,yd,A,B,OCP)
save('u.mat','u')

end
