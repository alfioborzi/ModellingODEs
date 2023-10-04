%Lorenz-Modell
% d/dt x = sigma(y-x)
% d/dt y = rx-y-xz
% d/dt z = xy - bz
% x = [x;y;z]
function xdot = lorenz(x,sigma,r,b)
	xdot = [0; 0; 0];
	xdot(1) = sigma*(x(2)-x(1));
	xdot(2) = r*x(1)-x(2)-x(1)*x(3);
	xdot(3) = x(1)*x(2)-b*x(3);
end
