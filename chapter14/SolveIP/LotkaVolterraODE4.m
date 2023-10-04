% the Lotka-Volterra model
%
function xdot = LotkaVolterraODE4(x,a,b,c,d)
	
	xdot = [x(1) * ( a - b * x(2) ) ;
			x(2) * ( c * x(1) - d )];
end
