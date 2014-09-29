%------------ MODEL --------------%
function y = AckermannModel(x,u,dT,L)
y(1,1) = x(1) + dT*u(1)*cos(x(3));
y(2,1) = x(2) + dT*u(1)*sin(x(3));
y(3,1) = x(3) + dT*u(1)/L*tan(u(2));

% y(1,1) = x(1) + u(1)*cos(x(3));
% y(2,1) = x(2) + u(1)*sin(x(3));
% y(3,1) = x(3) + u(1)/L*tan(u(2));
