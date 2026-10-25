clc
clear all
format long
%% Inicialização Código
%Parâmetros gerais
global mu                   % Definindo variável global
mu = 3.986004418e14;        % Constante
Earth_r = 6378140;          % Raio da Terra
a0 = Earth_r + 600e3;       % Distância do Líder
n0 = sqrt(mu/a0^3);         % 
e0 = 0.1;                   % Ecêntricidade
omega0 = 0;                 % 
f0 = 0;                     % 
theta = omega0 + f0;        % Argumento de latitude
r = 10*sqrt(2);             % Ruído das medidas
passo = 5;                  % Passo    
num_orbits = 4;
T = num_orbits*5752;        % Período
h = T/passo;                % Pontos de amostragem

%% Condições iniciais
xn = -0.01127;
yn = 0;
zn = 0.1;
xdot_n = 0.02;
ydot_n = 0.02;
zdot_n = 0;
r0 = (a0*(1-e0^2))/(1+e0*cos(f0));
r0dot = e0*sin(f0)*sqrt(mu/(a0*(1-e0^2)));
thetadot = (sqrt(mu/(a0*(1-e0^2))^3))*(1+e0*cos(f0))^2;

% Valores iniciais para vetor de estados
x1o = xn*a0            % x                 [m]
x2o = yn*a0            % y                 [m]
x3o = zn*a0            % z                 [m]
x4o = xdot_n*(n0*a0)   % xdot              [m/s]
x5o = ydot_n*(n0*a0)   % ydot              [m/s]
x6o = zdot_n*(n0*a0)   % zdot              [m/s]
x7o = r0               % r0                [m]
x8o = r0dot            % r0dot             [m]
x9o = theta            % theta0            [rad]
x10o = thetadot        % theta0dot         [rad/s]

% Valores iniciais para matriz de covariancia
p1o = 20;               % sigma_x            [m]
p2o = 20;               % sigma_y            [m]
p3o = 20;               % sigma_z            [m]
p4o = 0.1;              % sigma_vx           [m/s]
p5o = 0.1;              % sigma_vy           [m/s]
p6o = 0.1;              % sigma_vz           [m/s]
p7o = 1.5e-1;           % sigma_r0           [m]
p8o = 1e-10;            % sigma_r0dot        [m/s]
p9o = 1e-10*pi/180;     % sigma_theta0       [rad]
p10o = 0.001*pi/180;    % sigma_thata0dot    [rad/s]   

% Vetor de estado e matriz de covariância
Xo = [x1o; x2o; x3o; x4o; x5o; x6o; x7o; x8o; x9o; x10o]; % [m]
Po = diag([p1o^2; p2o^2; p3o^2; p4o^2; p5o^2; p6o^2; ...
    p7o^2; p8o^2; p9o^2; p10o^2]);                        % [m]

Xc = Xo;
Pc = Po;

% Condiçõe iniciais
CI = [Xc; Pc(:,1); Pc(:,2); Pc(:,3); Pc(:,4); Pc(:,5); Pc(:,6); Pc(:,7); Pc(:,8); Pc(:,9); Pc(:,10)];

tab1 = zeros(floor(h),111);

for i = 1:h,
    %Propagação da mecânica
    [t,Xp] = ode45(@edo,[(passo*i-passo) i*passo],CI);
    CI = Xp(end,:);
    tab1(i,:) = [i*passo CI];
end

Xp = tab1;
S = Xp(:,2:4);
d = size(Xp);
R = (r*randn(d(1),3));

data = [Xp(:,1) (S+R)];
data_verdade = [Xp(:,1) (Xp(:,2:11))];
dados_hiperdim = [Xp(:,1) (S+R*5e2)];

save data.dat data -ascii
save data_verdade.dat data_verdade -ascii
save dados_hiperdim.dat dados_hiperdim -ascii