function Xp = edo(t,ci)

global mu

q1o = 1;                  % q_x             [m]
q2o = 1;                  % q_y             [m]
q3o = 1;                  % q_z             [m]
q4o = 0.01;               % q_vx            [m/s]
q5o = 0.01;               % q_vy            [m/s]
q6o = 0.01;               % q_vz            [m/s]
q7o = 10e-6;                  % r0              [m]
% q8o = 0.01;              % r0dot           [m/s]
% q9o = 1e-5*pi/180;       % theta0          [rad]
% q10o = 1e-4*pi/180;       % theta0dot       [rad/s]   
q8o = 1e-20;              % r0dot           [m/s]
q9o = 1e-20*pi/180;       % theta0          [rad]
q10o = 1e-20*pi/180;       % theta0dot       [rad/s]   

Q = diag([q1o^2; q2o^2; q3o^2; q4o^2; q5o^2; q6o^2; q7o^2; q8o^2; q9o^2; q10o^2]);

x1 = ci(1);               % x               [m]
x2 = ci(2);               % y               [m]
x3 = ci(3);               % z               [m]
x4 = ci(4);               % xdot            [m/s]
x5 = ci(5);               % ydot            [m/s]
x6 = ci(6);               % zdot            [m/s]
x7 = ci(7);               % r0              [m]
x8 = ci(8);               % r0dot           [m/s]
x9 = ci(9);               % theta           [rad]
x10 = ci(10);             % thetadot        [rad/s]
    
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

F = [0,0,0,1,0,0,0,0,0,0;
     0,0,0,0,1,0,0,0,0,0;
     0,0,0,0,0,1,0,0,0,0;
     theta0dot^2-mu/((r0+x)^2+y^2+z^2)^(3/2)+(3*mu*(2*r0+2*x)*(r0+x))/(2*((r0+x)^2+y^2+z^2)^(5/2)),...
     (3*mu*y*(r0+x))/((r0+x)^2+y^2+z^2)^(5/2)-(2*r0dot*theta0dot)/r0,...
     (3*mu*z*(r0+x))/((r0+x)^2+y^2+z^2)^(5/2),0,2*theta0dot,0,...
     (3*mu*(2*r0+2*x)*(r0+x))/(2*((r0+x)^2+y^2+z^2)^(5/2))-mu/((r0+x)^2+y^2+z^2)^(3/2)-(2*mu)/r0^3+(2*r0dot*theta0dot*y)/r0^2,...
     -(2*theta0dot*y)/r0,0,2*ydot+2*theta0dot*x-(2*r0dot*y)/r0;
     (2*r0dot*theta0dot)/r0+(3*mu*y*(2*r0+2*x))/(2*((r0+x)^2+y^2+z^2)^(5/2)),...
     theta0dot^2-mu/((r0+x)^2+y^2+z^2)^(3/2)+(3*mu*y^2)/((r0+x)^2+y^2+z^2)^(5/2),(3*mu*y*z)/((r0+x)^2+y^2+z^2)^(5/2),...
     -2*theta0dot,0,0,(3*mu*y*(2*r0+2*x))/(2*((r0+x)^2+y^2+z^2)^(5/2))-(2*r0dot*theta0dot*x)/r0^2,...
     (2*theta0dot*x)/r0, 0, 2*theta0dot*y-2*xdot+(2*r0dot*x)/r0;
     (3*mu*z*(2*r0+2*x))/(2*((r0+x)^2+y^2+z^2)^(5/2)),(3*mu*y*z)/((r0+x)^2+y^2+z^2)^(5/2),...
     (3*mu*z^2)/((r0+x)^2+y^2+z^2)^(5/2)-mu/((r0+x)^2+y^2+z^2)^(3/2),0,0,0,(3*mu*z*(2*r0+2*x))/(2*((r0+x)^2+y^2+z^2)^(5/2)),0,0,0;
     0,0,0,0,0,0,0,1,0,0;
     0,0,0,0,0,0,(2*mu)/r0^3+theta0dot^2,0,0,2*r0*theta0dot;
     0,0,0,0,0,0,0,0,0,1;
     0,0,0,0,0,0,(2*r0dot*theta0dot)/r0^2,-(2*theta0dot)/r0,0,-(2*r0dot)/r0];


dx7 = r0dot;                                                %r0dot
dx8 = r0*theta0dot^2-mu/(r0^2);                             %r0dotdot
dx9 = theta0dot;                                            %theta0dot
dx10 = -2*r0dot*theta0dot/r0;                               %theta0dotdot
dx1 = xdot;                                                 %xdot
dx2 = ydot;                                                 %ydot
dx3 = zdot;                                                 %zdot 
dx4 = 2*theta0dot*ydot+dx10*y+(theta0dot^2)*x...
      -mu*(r0+x)/(((r0+x)^2+y^2+z^2)^(3/2))+mu/(r0^2);      %xdotdot
dx5 = -2*theta0dot*xdot-dx10*x+(theta0dot^2)*y ...
      -mu*y/(((r0+x)^2+y^2+z^2)^(3/2));                     %ydotdot
dx6 = -mu*z/(((r0+x)^2+y^2+z^2)^(3/2));                     %zdotdot

Pp = F*p+p*F'+Q;

dx = [dx1; dx2; dx3; dx4; dx5; dx6; dx7; dx8; dx9; dx10];
dp = [Pp(:,1); Pp(:,2); Pp(:,3); Pp(:,4); Pp(:,5); Pp(:,6); Pp(:,7); Pp(:,8); Pp(:,9); Pp(:,10)];
Xp = [dx; dp];