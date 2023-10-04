function J = modelTG(b)
% computation of the error functional of the ODE model
% b - vector of parameters

global xi T N tdata xdata x0

%% ODE model 
% 
    function dx = f(t,x)
        dx = zeros(2,1);      

dx(1) = -xi*x(1)*log(x(1)/x(2));
dx(2) = (b(1)*x(1)-b(2)*x(1)^(2/3)*x(2));
    end

%% numerical integration set up
N = 1000; 
dt = T/N; 

tspan = [0:dt:max(tdata)];
[tsol,xsol] = ode45(@f,tspan,x0);

%% plot results

ltype = {'b-','r--','m-.'}; 

figure(1)
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
drawnow

%% Predicted values x(tdata)

xpred = interp1(tsol,xsol,tdata);

%% Error functional

J = 0;
for i = 1:length(tdata)
    J = J + sum((xpred(i,:)-xdata(i,:)).^2);
end

end