clc
clear all

global mu
format long
r = 10; % erro da medida
passo = 5; % amostragem das medidas (s)

% T = 2*pi*sqrt(a0^3/mu);
T = 4*5752;

h = T/passo;


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

%Valores iniciais para matriz de covariancia
p1o = 10;
p2o = 10; 
p3o = 10; 
p4o = 0.1;
p5o = 0.1; 
p6o = 0.1; 
p7o = 0; 
p8o = 0; 
p9o = 0;
p10o = 0;  

Xo = [x1o; x2o; x3o; x4o; x5o; x6o; x7o; x8o; x9o; x10o];
Po = diag([p1o^2; p2o^2; p3o^2; p4o^2; p5o^2; p6o^2; p7o^2; p8o^2; p9o^2; p10o^2]);

Xc = Xo;
Pc = Po;

CI = [Xc; Pc(:,1); Pc(:,2); Pc(:,3); Pc(:,4); Pc(:,5); Pc(:,6); Pc(:,7); Pc(:,8); Pc(:,9); Pc(:,10)];

tab1 = zeros(floor(h),111);

% for i = 1:h(1),
for i = 1:h,
    %Fase Propagação
    [t,Xp] = ode45(@edo,[(passo*i-passo) i*passo],CI);
    CI = Xp(end,:);
    tab1(i,:) = [i*passo CI];
    
end
% [t,Xp] = ode45(@edo,0:passo:(5752*4),CI);

Xp = tab1;

S = Xp(:,2:4);
d = size(Xp);
R = (r*randn(d(1),3));

data = [Xp(:,1) (S+R)];


save data.dat data -ascii

load data.dat
xyz = data(:,2:4);



figure(1)
subplot(2,2,1)
plot(xyz(:,1)/a0,xyz(:,2)/a0,'LineWidth',1)
hold on
plot(S(:,1)/a0,S(:,2)/a0)
xlabel('x (m)');
ylabel('y (m)');
xlim([-0.06 0.06]);
ylim([-0.08 0.08]);
legend('Órbita Medida','Órbita Estimada','Location','NorthWest');
grid on
subplot(2,2,2)
plot(xyz(:,1)/a0,xyz(:,3)/a0,'LineWidth',1)
hold on
plot(S(:,1)/a0,S(:,3)/a0)
xlabel('x (m)');
ylabel('z (m)');
xlim([-0.12 0.12]);
ylim([-0.14 0.14]);
grid on
subplot(2,2,3)
plot(xyz(:,2)/a0,xyz(:,3)/a0,'LineWidth',1)
hold on
plot(S(:,2)/a0,S(:,3)/a0)
xlabel('y (m)');
ylabel('z (m)');
xlim([-0.12 0.12]);
ylim([-0.14 0.14]);
grid on
subplot(2,2,4)
plot3(xyz(:,1)/a0,xyz(:,2)/a0,xyz(:,3)/a0)
hold on
plot3(S(:,1)/a0,S(:,2)/a0,S(:,3)/a0,'LineWidth',1)
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
xlim([-0.12 0.12]);
ylim([-0.12 0.12]);
zlim([-0.12 0.12]);
hold on
plot3(xyz(1,1)/a0,xyz(1,2)/a0,xyz(1,3)/a0,'*r','LineWidth',1)
grid on
