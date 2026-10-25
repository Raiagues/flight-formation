%clc
%clear all

global mu
format long
r = 10^2; % erro da medida
passo = 5; % amostragem das medidas (s)

Earth_r = 6378140;
mu = 3.986004418e14;

a0 = Earth_r + 600e3;
e0 = 0.1;
w0 = 0;
f0 = 0;
teta0 = w0 + f0;

num_orbits = 4;

xn = -0.01127;
yn = 0;
zn = 0.1;
xdot_n = 0.02;
ydot_n = 0.02;
zdot_n = 0;

n0 = sqrt(mu/a0^3);

x_abs = xn*a0;
y_abs = yn*a0;
z_abs = zn*a0;
xdot_abs = xdot_n*(n0*a0);
ydot_abs = ydot_n*(n0*a0);
zdot_abs = zdot_n*(n0*a0);

T = 2*pi*sqrt(a0^3/mu);

r0 = (a0*(1-e0^2))/(1+e0*cos(f0));
r0dot = e0*sin(f0)*sqrt(mu/(a0*(1-e0^2)));
thetadot = (sqrt(mu/(a0*(1-e0^2))^3))*(1+e0*cos(f0))^2;

x1o = x_abs;
x2o = y_abs;
x3o = z_abs;
x4o = xdot_abs;
x5o = ydot_abs;
x6o = zdot_abs;
x7o = r0;
x8o = r0dot;
x9o = teta0;
x10o = thetadot;

p1o = 1;
p2o = 1; 
p3o = 1; 
p4o = 0.001;
p5o = 0.001; 
p6o = 0.001; 
p7o = 0; 
p8o = 0; 
p9o = 0;
p10o = 0;

Xo = [x1o; x2o; x3o; x4o; x5o; x6o; x7o; x8o; x9o; x10o];
Po = diag([p1o^2; p2o^2; p3o^2; p4o^2; p5o^2; p6o^2; p7o^2; p8o^2; p9o^2; p10o^2]);

Xc = Xo;
Pc = Po;

CI = [Xc; Pc(:,1); Pc(:,2); Pc(:,3); Pc(:,4); Pc(:,5); Pc(:,6); Pc(:,7); Pc(:,8); Pc(:,9); Pc(:,10)];

[t,Xp] = ode45(@edo,0:passo:(5752*2),CI);

S = Xp(:,1:3);
d = size(Xp);
R = (r*randn(d(1),3));

data = [t (S+R)];

% tf = round(T*num_orbits);
% 
% s0 = [x_abs; y_abs; z_abs; xdot_abs; ydot_abs; zdot_abs; r0; r0dot; teta0; thetadot];
% 
% [t,s] = ode45(@nonlinear_rme_unp,[0 tf],s0);

%S = s;

%r = 100^2;
% r = 0;
S = Xp(:,1:3);
d = size(Xp);
R = (r*randn(d(1),3));

data = [t (S+R)];

save data.dat data -ascii
% 
% load data.dat
% xyz = data(:,2:4);


figure(1)
subplot(2,2,1)
plot(xyz(:,1),xyz(:,2),'LineWidth',1)
hold on
plot(S(:,1),S(:,2))
xlabel('x (m)');
ylabel('y (m)');
legend('Órbita Medida','Órbita Estimada','Location','NorthWest');
grid on
subplot(2,2,2)
plot(xyz(:,1),xyz(:,3),'LineWidth',1)
hold on
plot(S(:,1),S(:,3))
xlabel('x (m)');
ylabel('z (m)');
grid on
subplot(2,2,3)
plot(xyz(:,2),xyz(:,3),'LineWidth',1)
hold on
plot(S(:,2),S(:,3))
xlabel('y (m)');
ylabel('z (m)');
grid on
subplot(2,2,4)
plot3(xyz(:,1),xyz(:,2),xyz(:,3))
hold on
plot3(S(:,1),S(:,2),S(:,3),'LineWidth',1)
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
hold on
plot3(xyz(1,1),xyz(1,2),xyz(1,3),'*r','LineWidth',1)
grid on
