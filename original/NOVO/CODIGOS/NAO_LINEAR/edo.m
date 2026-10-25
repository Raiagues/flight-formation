function Xp = edo(t,ci)

global mu

w1 = 0;
w2 = 0;
w3 = 0;

G = [0 0 0;0 0 0;0 0 0;1 0 0;0 1 0;0 0 1;0 0 0;0 0 0;0 0 0;0 0 0];
Q = [w1^2 0 0;0 w2^2 0;0 0 w3^2];

%x1 = ci(1); % x
%x2 = ci(2); % y
%x3 = ci(3); % z
%x4 = ci(4); % xdot
%x5 = ci(5); % ydot
%x6 = ci(6); % zdot
%x7 = ci(7); % r0
%x8 = ci(8); % r0dot
%x9 = ci(9); % theta
%x10 = ci(10); % thetadot

x1 = ci(1); % x
x2 = ci(2); % y
x3 = ci(3); % z
x4 = ci(4); % xdot
x5 = ci(5); % ydot
x6 = ci(6); % zdot
x7 = ci(7); % r0
x8 = ci(8); % r0dot
x9 = ci(9); % theta
x10 = ci(10); % thetadot

    
p = [ci(11:20) ci(21:30) ci(31:40) ci(41:50) ci(51:60) ci(61:70) ci(71:80) ci(81:90) ci(91:100) ci(101:110)];

x = x1;
y = x2;
z = x3;
xdot = x4;
ydot = x5;
zdot = x6;
r0 = x7;
r0dot = x8;
theta0dot = x10;

F = [0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
     theta0dot^2 - mu/((r0 + x)^2 + y^2 + z^2)^(3/2) + (3*mu*(2*r0 + 2*x)*(r0 + x))/(2*((r0 + x)^2 + y^2 + z^2)^(5/2)), (3*mu*y*(r0 + x))/((r0 + x)^2 + y^2 + z^2)^(5/2) - (2*r0dot*theta0dot)/r0, (3*mu*z*(r0 + x))/((r0 + x)^2 + y^2 + z^2)^(5/2), 0, 2*theta0dot, 0, (3*mu*(2*r0 + 2*x)*(r0 + x))/(2*((r0 + x)^2 + y^2 + z^2)^(5/2)) - mu/((r0 + x)^2 + y^2 + z^2)^(3/2) - (2*mu)/r0^3 + (2*r0dot*theta0dot*y)/r0^2, -(2*theta0dot*y)/r0, 0, 2*ydot + 2*theta0dot*x - (2*r0dot*y)/r0;
     (2*r0dot*theta0dot)/r0 + (3*mu*y*(2*r0 + 2*x))/(2*((r0 + x)^2 + y^2 + z^2)^(5/2)), theta0dot^2 - mu/((r0 + x)^2 + y^2 + z^2)^(3/2) + (3*mu*y^2)/((r0 + x)^2 + y^2 + z^2)^(5/2), (3*mu*y*z)/((r0 + x)^2 + y^2 + z^2)^(5/2), -2*theta0dot, 0, 0, (3*mu*y*(2*r0 + 2*x))/(2*((r0 + x)^2 + y^2 + z^2)^(5/2)) - (2*r0dot*theta0dot*x)/r0^2, (2*theta0dot*x)/r0, 0, 2*theta0dot*y - 2*xdot + (2*r0dot*x)/r0;
     (3*mu*z*(2*r0 + 2*x))/(2*((r0 + x)^2 + y^2 + z^2)^(5/2)), (3*mu*y*z)/((r0 + x)^2 + y^2 + z^2)^(5/2), (3*mu*z^2)/((r0 + x)^2 + y^2 + z^2)^(5/2) - mu/((r0 + x)^2 + y^2 + z^2)^(3/2), 0, 0, 0, (3*mu*z*(2*r0 + 2*x))/(2*((r0 + x)^2 + y^2 + z^2)^(5/2)), 0, 0, 0;
     0, 0, 0, 0, 0, 0, 0, 1, 0, 0;
     0, 0, 0, 0, 0, 0, (2*mu)/r0^3 + theta0dot^2, 0, 0, 2*r0*theta0dot;
     0, 0, 0, 0, 0, 0, 0, 0, 0, 1;
     0, 0, 0, 0, 0, 0, (2*r0dot*theta0dot)/r0^2, -(2*theta0dot)/r0, 0, -(2*r0dot)/r0];


dx7 = r0dot;                   %r0dot
dx8 = r0*theta0dot^2 - mu/(r0^2);   %r0dotdot
dx9 = theta0dot;                  %theta0dot
dx10 = -2*r0dot*theta0dot/r0;        %theta0dotdot

dx1 = xdot;
dx2 = ydot;
dx3 = zdot;
dx4 = 2*theta0dot*ydot + dx10*y + (theta0dot^2)*x - mu*(r0+x)/((r0+x)^2+y^2+z^2)^(3/2) + mu/(r0^2);
dx5 = -2*theta0dot*xdot - dx10*x + (theta0dot^2)*y - mu*y/((r0+x)^2+y^2+z^2)^(3/2);
dx6 = -mu*z/((r0+x)^2+y^2+z^2)^(3/2);

Pp = F*p + p*F' + G*Q*G';

dx = [dx1; dx2; dx3; dx4; dx5; dx6; dx7; dx8; dx9; dx10];
dp = [Pp(:,1); Pp(:,2); Pp(:,3); Pp(:,4); Pp(:,5); Pp(:,6); Pp(:,7); Pp(:,8); Pp(:,9); Pp(:,10)];
Xp = [dx; dp];