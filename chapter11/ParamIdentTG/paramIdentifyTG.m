function paramIdentifyTG

% Parameters identification of a tumor growth ODE model to data
% the model and the error functional are defined in modelTG.m

clearvars -global
global xi T N tdata xdata x0

%% Time horizon

T = 30; 

%% initial conditions p0 and q0

x0(1) = 80;
x0(2) = 100;

%% parameter values for generating data

xi=0.084;

% Values of the parameters b = b(1) and d = b(2)
b(1)=5.85;
b(2)=0.00873;

%% generate data using the model

    function dx = fg(t,x)
        dx = zeros(2,1);
dx(1) = -xi*x(1)*log(x(1)/x(2));
dx(2) = (b(1)*x(1)-b(2)*x(1)^(2/3)*x(2));
    end

%% numerical integration 

N = 1000; 
dt = T/N; 
tspan = [0:dt:T];
[tsol,xsol] = ode45(@fg,tspan,x0);

%% Number of measurements : kd
[nd,md] = size(tsol);
kd = 5;
sd = round(nd/kd);
iData = 1:sd:nd; 

%% add noise

for id=1:kd
    tdata(id)   = tsol(iData(id)); 
    xdata(id,1) = xsol(iData(id),1) * (1 + 0.1 * (rand - 0.5)) ;  
    xdata(id,2) = xsol(iData(id),2) * (1 + 0.1 * (rand - 0.5)) ;
end


%% initial guess of parameter values

b(1)=0.5;
b(2)=0.5;

%% minimization step
%% use fmincon to add bounds or set J=infinity in modelTG if b out of bounds

 [bmin, Smin] = fminsearch(@modelTG,b);


%%
ltype = {'b-','r--','m-.'}; 


figure(2)
    subplot(1,2,1)
    plot(tdata,xdata(:,1),'o','MarkerSize',10);
    hold on
    plot(tsol,xsol(:,1),ltype{1},'Linewidth',2);
    hold off
    ylabel(['p']); xlabel(['days']); 

    subplot(1,2,2)
    plot(tdata,xdata(:,2),'o','MarkerSize',10);
    hold on
    plot(tsol,xsol(:,2),ltype{1},'Linewidth',2);
    hold off
    ylabel(['q']); xlabel(['days']); 

print('-depsc2', 'estimationTG01.eps','-b0'); 
print('-dpdf', 'estimationTG01.pdf','-b0');


disp('Estimated parameters b and d:');
disp(bmin)
disp('Value of the error functional:');
disp(Smin)

end
