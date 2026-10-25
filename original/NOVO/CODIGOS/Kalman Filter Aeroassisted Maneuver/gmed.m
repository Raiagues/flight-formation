% Gerador de medidas
%
format long
r = 5e-4; % erro da medida
passo = 5; % amostragem das medidas (s)
%
% Condição inicial
% x1o = 6611.3586; % (X) km
% x2o = 0; % (Y) km
% x3o = 0; % (Z) km
% x4o = 0; % (Vx) km/s
% x5o = 7.803406; % (Vy) km/s
% x6o = 0; % (Vz) km/s
% x7o = 398600.4418; % (mi) km^3/s^2
% x8o = 2; % (Cd) adm
%
x1o = 6.498139999999999e3; % (X) km
x2o = 0; % (Y) km
x3o = 0; % (Z) km
x4o = 0; % (Vx) km/s
x5o = 8.076552708163819; % (Vy) km/s
x6o = 0; % (Vz) km/s
x7o = 398600.4418; % (mi) km^3/s^2
x8o = 2; % (Cd) adme3; % (X) km
%
p1o = 1; % (Sigma X) km
p2o = 1; % (Sigma Y) km
p3o = 1; % (Sigma Z) km
p4o = 0.001; % (Sigma Vx) km/s
p5o = 0.001; % (Sigma Vy) km/s
p6o = 0.001; % (Sigma Vz) km/s
p7o = 0.1; % (Sigma mi) km^3/s^2
p8o = 0.2; % (Sigma Cd) adm
%
Xo = [x1o; x2o; x3o; x4o; x5o; x6o; x7o; x8o];
Po = diag([p1o^2; p2o^2; p3o^2; p4o^2; p5o^2; p6o^2; p7o^2; p8o^2]);
%
Xc = Xo;
Pc = Po;
CI = [Xc; Pc(:,1); Pc(:,2); Pc(:,3); Pc(:,4); Pc(:,5); Pc(:,6);...
    Pc(:,7); Pc(:,8)];
%
[t,Xp] = ode45(@propag,0:passo:(5752*2),CI);
%
%h = sqrt(Xp(:,1).^2 + Xp(:,2).^2 + Xp(:,3).^2) - 6378.14;
%
S = Xp(:,1:3);
%
d = size(Xp);
R = (r*randn(d(1),3));
%
data = [t (S + R)];
%
save data.dat data -ascii

%
% Nota: Matriz data = 1087 linhas x 4 colunas