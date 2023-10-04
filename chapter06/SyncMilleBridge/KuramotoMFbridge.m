function der = KuramotoMFbridge(wtheta,omega,omega0,a,b,g,C,alpha)
%
% f = omega + Ac * sin(Psi - theta + alpha);

x = wtheta(1);
v = wtheta(2);

theta = wtheta(3:end);

Ac = C * sqrt( x^2 + v^2 / a );

if abs(v) > 1.0e-10
Psi = atan(omega0 * x / v );
else
Psi = pi/2;
end

force =  g * sum( sin( theta ) ) ;

        dxdt  = v;                      % set dx/dt = velocity
        dvdt  = -a * x - b * v + force;  % set dv/dt = acceleration
        dthdt = omega + Ac * sin(Psi - theta + alpha);
        der = [dxdt; dvdt; dthdt];  % return the derivatives
 
end