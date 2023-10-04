function [J] = get_J(y,u,yd,OCP)
nu=OCP.nu;
dt=OCP.dt;
beta=OCP.beta;
J=0.5*((y-yd)'*(y-yd)) + (nu/2)*u*(u')*dt + beta*sum(abs(u))*dt; % erster Term y-y_d nur Endzeit
end

