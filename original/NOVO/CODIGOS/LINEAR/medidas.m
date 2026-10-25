clc
clear all
format long
%% Inicialização Código
%Parâmetros gerais
global mu                           % Definindo variável global
mu = 3.986004418e14;                % Parametro gravitacional da terra
Earth_r = 6378140;                  % Raio da Terra
a0 = Earth_r + 600e3;               % Distância do Líder
r = 10*sqrt(2);                     % Ruído das medidas
passo = 5;                          % Passo    
num_orbits = 4;                     % Numero de orbitas
T = num_orbits*5752;  % Período
h = T/passo;                        % Pontos de amostragem

%% Condições iniciais
x1o = 69.78e3;          % x                 [m]
x2o = 139.56e3;         % y                 [m]
x3o = 104.67e3;         % z                 [m]
x4o = 7.5579;           % xdot              [m/s]
x5o = -151.116;         % ydot              [m/s]
x6o = 15.116;           % zdot              [m/s]

% Valores iniciais para matriz de covariancia
p1o = 20;               % sigma_x           [m]
p2o = 20;               % sigma_y           [m]
p3o = 20;               % sigma_z           [m]
p4o = 0.1;              % sigma_vx          [m/s]
p5o = 0.1;              % sigma_vy          [m/s]
p6o = 0.1;              % sigma_vz          [m/s]

Xo = [x1o; x2o; x3o; x4o; x5o; x6o];
Po = diag([p1o^2; p2o^2; p3o^2; p4o^2; p5o^2; p6o^2]);

Xc = Xo;
Pc = Po;
CI = [Xc; Pc(:,1); Pc(:,2); Pc(:,3); Pc(:,4); Pc(:,5); Pc(:,6)];

tab1 = zeros(floor(h),43);

for i = 1:h,
    %Propagação da mecanica
    [t,Xp] = ode45(@edo,[(passo*i-passo) i*passo],CI);
    CI = Xp(end,:);
    tab1(i,:) = [i*passo CI];
end


Xp = tab1;
S = Xp(:,2:4);
d = size(Xp);
R = (r*randn(d(1),3));

data = [Xp(:,1) (S+R)];
data_verdade = [Xp(:,1) (Xp(:,2:7))];
dados_hiperdim = [Xp(:,1) (S+R*100)];

save data.dat data -ascii
save data_verdade.dat data_verdade -ascii
save dados_hiperdim.dat dados_hiperdim -ascii

% [t,Xp] = ode45(@edo,0:passo:(5752*2),CI);
% 
% S = Xp(:,1:3);
% d = size(Xp);
% R = (r*randn(d(1),3));
% 
% data = [t (S+R)];

%data = [t (S)];

% save data.dat data -ascii
%
% Nota: Matriz data = 1087 linhas x 4 colunas