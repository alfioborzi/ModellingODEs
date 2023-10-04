% Modified Prey Predator Model
% based on the Lotka-Volterra model
%
function xdot = preypredatorODE(t,x,a,b,d)
	
	xdot = [x(1) * (  (1 - x(1)) -  a * x(2)/(x(1) + d) );
			x(2) * b * ( 1 - x(2)/x(1) )];
end
