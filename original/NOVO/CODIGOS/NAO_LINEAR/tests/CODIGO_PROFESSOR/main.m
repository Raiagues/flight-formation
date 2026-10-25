clc
close all
clear all
%% Inicialização Código
% Parâmetros gerais
global mu                   % Definindo variável global
mu = 3.986004418e14;        % Constante
Earth_r = 6378140;          % Raio da Terra
a0 = Earth_r + 600e3;       % Distância do Líder
n0 = sqrt(mu/a0^3);         %
e0 = 0.1;                   % Ecêntricidade
omega0 = 0;                 % 
f0 = 0;                     %
theta = omega0 + f0;        % Argumento de latitude

% Importando medidas iniciais
format long
load data.dat               % Coordenadas medidas
load data_verdade.dat       % Coordenadas sem ruídos
load dados_hiperdim.dat     % Coordenadas com ruídos hiperdimensionados

% Declarando coordenadas sem ruídos
xyz_verdade = data_verdade(:,2:11);         % Coordenadas sem ruídos
xyz = data(:,2:4);                          % Coordenadas medidas
xyz_hiperdim = dados_hiperdim(:,2:4);       % Coordenadas com ruídos hiperdimensionados
passo  = 5;                                 % Passo de amostragem
h = size(xyz);

%Definindo matrizes para guardar resultados
tabt = zeros(h(1),1);                       % Tempo de simulação
tab1 = zeros(h(1),10);                      % Estado estimado
tab2 = zeros(h(1),6);                       % Sequência de inovação
tab3 = zeros(h(1),10);                      % Desvios
tab4 = zeros(h(1),10);                      % Estado propagado

%% Condições iniciais
xn = -0.01127;                                              % [adm]
yn = 0;                                                     % [adm]
zn = 0.1;                                                   % [adm]
xdot_n = 0.02;                                              % [adm]
ydot_n = 0.02;                                              % [adm]
zdot_n = 0;                                                 % [adm]
r0 = (a0*(1-e0^2))/(1+e0*cos(f0));                          % [m]
r0dot = e0*sin(f0)*sqrt(mu/(a0*(1-e0^2)));                  % [m/s]
thetadot = (sqrt(mu/(a0*(1-e0^2))^3))*(1+e0*cos(f0))^2;     % [rad]

% Valores iniciais para vetor de estados
x1o = xn*a0;            % x                 [m]
x2o = yn*a0;            % y                 [m]
x3o = zn*a0;            % z                 [m]
x4o = xdot_n*(n0*a0);   % xdot              [m/s]
x5o = ydot_n*(n0*a0);   % ydot              [m/s]
x6o = zdot_n*(n0*a0);   % zdot              [m/s]
x7o = r0;               % r0                [m]
x8o = r0dot;            % r0dot             [m]
x9o = theta;            % theta0            [rad]
x10o = thetadot;        % theta0dot         [rad/s]

% Valores iniciais para matriz de covariancia
p1o = 20;               % sigma_x            [m]
p2o = 20;               % sigma_y            [m]
p3o = 20;               % sigma_z            [m]
p4o = 0.1;              % sigma_vx           [m/s]
p5o = 0.1;              % sigma_vy           [m/s]
p6o = 0.1;              % sigma_vz           [m/s]
% p7o = 1.5e-1;           % sigma_r0           [m]
p7o = 1e-20;           % sigma_r0           [m]
% p8o = 1e-10;            % sigma_r0dot        [m/s]
p8o = 1e-20;            % sigma_r0dot        [m/s]
% p9o = 1e-10*pi/180;     % sigma_theta0       [rad]
% p10o = 0.001*pi/180;    % sigma_thata0dot    [rad/s]   
p9o = 1e-20*pi/180;     % sigma_theta0       [rad]
p10o = 1e-20*pi/180;    % sigma_thata0dot    [rad/s]   

%% Inicialização do Filtro de Kalman
% Vetor de estado e matriz de covariância
Xo = [x1o; x2o; x3o; x4o; x5o; x6o; x7o; x8o; x9o; x10o]; % [m]
Po = diag([p1o^2; p2o^2; p3o^2; p4o^2; p5o^2; p6o^2; ...
    p7o^2; p8o^2; p9o^2; p10o^2]);                        % [m]

% Matriz erro de observação
R = diag([(10*sqrt(2))^2 (10*sqrt(2))^2 (10*sqrt(2))^2]); % [m]

% Matrizes condições iniciais
xb = Xo;                                                  % [m]
Pb = Po;                                                  % [m]

% Matriz de observação
H = [1 0 0 0 0 0 0 0 0 0;0 1 0 0 0 0 0 0 0 0;0 0 1 0 0 0 0 0 0 0];

%% Filtro de Kalman
for i = 1:h(1),
    
    % Fase Atualização
    K = Pb*H'*inv(H*Pb*H'+R);                       
    Pc = (eye(10,10)-K*H)*Pb;
    y = [xyz(i,1); xyz(i,2); xyz(i,3)];
    y_verdade = [xyz_verdade(i,1); xyz_verdade(i,2); xyz_verdade(i,3)];
    
    % Calculando sequência de inovação
    res = y-H*xb;
    Xc = xb+K*res;
    res2 = Xc(1:3)-xb(1:3);
    dvp = sqrt(diag(Pc));
    
    % Atualizando condições iniciais
    CI = [Xc; Pc(:,1); Pc(:,2); Pc(:,3); Pc(:,4); Pc(:,5); Pc(:,6); Pc(:,7); Pc(:,8); Pc(:,9); Pc(:,10)];
    
    % Fase Propagação
    [t,Xp] = ode45(@edo,[(passo*i-passo) i*passo],CI);
    d = size(Xp);
    Xb = Xp(d(1),:)';
    xb = Xb(1:10);
    Pb = [Xb(11:20) Xb(21:30) Xb(31:40) Xb(41:50) Xb(51:60) Xb(61:70) Xb(71:80) Xb(81:90) Xb(91:100) Xb(101:110)];
    
    % Guardando variáveis
    tabt(i,:) = i*passo;
    tab1(i,:) = [Xc(1),Xc(2),Xc(3),Xc(4),Xc(5),Xc(6),Xc(7),Xc(8),Xc(9),Xc(10)];
    tab2(i,:) = [res',res2'];
    tab3(i,:) = [dvp(1),dvp(2),dvp(3),dvp(4),dvp(5),dvp(6),dvp(7),dvp(8),dvp(9),dvp(10)];
    tab4(i,:) = [xb(1),xb(2),xb(3),xb(4),xb(5),xb(6),xb(7),xb(8),xb(9),xb(10)];
end

%% Organizando dados para plotagem
tm = cat(1,0,tabt);
tabx0 = [x1o x2o x3o x4o x5o x6o x7o x8o x9o x10o];
tab1 = cat(1,tabx0,tab1);
tabdv = tab3;
tabp0 = [p1o p2o p3o p4o p5o p6o p7o p8o p9o p10o];
tab3 = cat(1,tabp0,tab3);
smax_pos = [abs(3*tabdv(:,1)) abs(3*tabdv(:,2)) abs(3*tabdv(:,3))];
smin_pos = -smax_pos;
tab4 = cat(2,smax_pos,smin_pos); 

% Retirando os ruídos do vetor de medidas
resx = (tab1(2:end,1)-xyz_verdade(:,1))*5e2;
resy = (tab1(2:end,2)-xyz_verdade(:,2))*5e2;
resz = (tab1(2:end,3)-xyz_verdade(:,3))*5e2;
resvx = (tab1(2:end,4)-xyz_verdade(:,4))*10;
resvy = (tab1(2:end,5)-xyz_verdade(:,5))*10;
resvz = (tab1(2:end,6)-xyz_verdade(:,6))*10;
resr0 = (tab1(2:end,7)-xyz_verdade(:,7));
resrdot0 = (tab1(2:end,8)-xyz_verdade(:,8));
restheta = (tab1(2:end,9)-xyz_verdade(:,9));
resthetadot = (tab1(2:end,10)-xyz_verdade(:,10));

% Criando um novo vetor com os ruídos hiperdimensinados
estado_aumentado = [xyz_verdade(:,1)+resx,xyz_verdade(:,2)+resy,xyz_verdade(:,3)+resz,...
                    xyz_verdade(:,4)+resvx,xyz_verdade(:,5)+resvy,xyz_verdade(:,6)+resvz,...
                    xyz_verdade(:,7)+resr0,xyz_verdade(:,8)+resrdot0,...
                    xyz_verdade(:,9)+restheta,xyz_verdade(:,10)+resthetadot];
                
%% Plotagens
% % Comparação entre órbita medida e órbita estimada
% figure(1)
% subplot(2,2,1)
% plot(xyz_hiperdim(2:end,1)/a0,xyz_hiperdim(2:end,2)/a0)
% hold on
% plot(estado_aumentado(2:end,1)/a0,estado_aumentado(2:end,2)/a0)
% xlim([-0.06 0.06]);
% ylim([-0.08 0.08]);
% xlabel('x/a0');
% ylabel('y/a0');
% legend('Órbita Medida','Órbita Estimada','Location','NorthWest');
% grid on
% subplot(2,2,2)
% plot(xyz_hiperdim(2:end,1)/a0,xyz_hiperdim(2:end,3)/a0)
% hold on
% plot(estado_aumentado(2:end,1)/a0,estado_aumentado(2:end,3)/a0)
% xlim([-0.12 0.12]);
% ylim([-0.14 0.14]);
% xlabel('x/a0');
% ylabel('z/a0');
% grid on
% subplot(2,2,3)
% plot(xyz_hiperdim(2:end,2)/a0,xyz_hiperdim(2:end,3)/a0)
% hold on
% plot(estado_aumentado(2:end,2)/a0,estado_aumentado(2:end,3)/a0)
% xlim([-0.12 0.12]);
% ylim([-0.14 0.14]);
% xlabel('y/a0');
% ylabel('z/a0');
% grid on
% subplot(2,2,4)
% plot3(xyz_hiperdim(2:end,1)/a0,xyz_hiperdim(2:end,2)/a0,xyz_hiperdim(2:end,3)/a0)
% hold on
% plot3(estado_aumentado(2:end,1)/a0,estado_aumentado(2:end,2)/a0,estado_aumentado(2:end,3)/a0)
% xlim([-0.12 0.12]);
% ylim([-0.12 0.12]);
% zlim([-0.12 0.12]);
% xlabel('x/a0');
% ylabel('y/a0');
% zlabel('z/a0');
% hold on
% plot3(xyz(1,1)/a0,xyz(1,2)/a0,xyz(1,3)/a0,'*r')
% grid on

% Desvio padrão da posição
figure(2)
desvio_x = tab3(:,1);
desvio_y = tab3(:,2);
desvio_z = tab3(:,3);
subplot(3,1,1)
semilogx(tm/3600, desvio_x, 'LineWidth',1);
grid on
xlabel('Tempo (h)');
ylabel('$\sigma_{x}$(m)', 'Interpreter', 'latex');
subplot(3,1,2)
semilogx(tm/3600, desvio_y, 'LineWidth',1);
grid on
xlabel('Tempo (h)');
ylabel('$\sigma_{y}$(m)', 'Interpreter', 'latex');
subplot(3,1,3)
semilogx(tm/3600, desvio_z, 'LineWidth',1);
grid on
xlabel('Tempo (h)');
ylabel('$\sigma_{z}$(m)', 'Interpreter', 'latex');

% Desvio padrão da velocidade
figure(3)
desvio_vx = tab3(:,4);
desvio_vy = tab3(:,5);
desvio_vz = tab3(:,6);
subplot(3,1,1)
semilogx(tm/3600, desvio_vx, 'LineWidth',1);
grid on
xlabel('Tempo (h)');
ylabel('$\sigma_{v_x}$(m/s)', 'Interpreter', 'latex');
subplot(3,1,2)
semilogx(tm/3600, desvio_vy, 'LineWidth',1);
grid on
xlabel('Tempo (h)');
ylabel('$\sigma_{v_y}$(m/s)', 'Interpreter', 'latex');
subplot(3,1,3)
semilogx(tm/3600, desvio_vz, 'LineWidth',1);
grid on
xlabel('Tempo (h)');
ylabel('$\sigma_{v_0}$(m/s)', 'Interpreter', 'latex');

% Desvio padrão da posição do líder e sua variação
figure(4)
desvio_r = tab3(:,7);
desvio_rdot = tab3(:,8);
subplot(2,1,1)
semilogx(tm/3600, desvio_r, 'LineWidth',1);
grid on
xlabel('Tempo (h)');
ylabel('$\sigma_{r_0}$(m)', 'Interpreter', 'latex');
subplot(2,1,2)
semilogx(tm/3600, desvio_rdot, 'LineWidth',1);
grid on
xlabel('Tempo (h)');
ylabel('$\sigma_{\dot{r}_0}$(m/s)', 'Interpreter', 'latex');

% Desvio padrão do argumento de latitude e sua variação
figure(5)
desvio_theta = tab3(:,9);
desvio_thetadot = tab3(:,10);
subplot(2,1,1)
semilogx(tm/3600, desvio_theta, 'LineWidth',1);
grid on
xlabel('Tempo (h)');
ylabel('$\sigma_{{\theta}_0}$(rad)', 'Interpreter', 'latex');
subplot(2,1,2)
semilogx(tm/3600, desvio_thetadot, 'LineWidth',1);
grid on
xlabel('Tempo (h)');
ylabel('$\sigma_{\dot{\theta}_0}$(rad/s)', 'Interpreter', 'latex');
% 
% % Comparação entre as posições nos instantes iniciais
% figure(6)
% subplot(3,1,1)
% plot(tm(1:end-2)/3600,xyz_hiperdim(1:end-1,1),'bs-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.4470 0.7410],...
%     'MarkerSize',3)
% hold on
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,1),'rd-','LineWidth',1,...
%     'MarkerFaceColor','r',...
%     'MarkerSize',3)
% hold on
% plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,1),'go-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.5 0],...
%     'MarkerSize',3)
% hold off
% set(findobj('color','g'),'Color',[0 0.5 0]);
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
% xlabel('Tempo (h)');
% ylabel('x (m)');
% xlim([0 120/3600])
% grid on
% legend('Medido','Estimado','Verdade', 'Location', 'southeast');
% lgd = legend;
% lgd.NumColumns = 3;
% subplot(3,1,2)
% plot(tm(1:end-2)/3600,xyz_hiperdim(1:end-1,2),'bs-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.4470 0.7410],...
%     'MarkerSize',3);
% hold on
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,2),'rd-','LineWidth',1,...
%     'MarkerFaceColor','r',...
%     'MarkerSize',3);
% hold on
% plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,2),'go-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.5 0],...
%     'MarkerSize',3);
% hold off
% set(findobj('color','g'),'Color',[0 0.5 0]);
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
% xlabel('Tempo (h)');
% ylabel('y (m)');
% xlim([0 120/3600])
% grid on
% subplot(3,1,3)
% plot(tm(1:end-2)/3600,xyz_hiperdim(1:end-1,3),'bs-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.4470 0.7410],...
%     'MarkerSize',3);
% hold on
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,3),'rd-','LineWidth',1,...
%     'MarkerFaceColor','r',...
%     'MarkerSize',3);
% hold on
% plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,3),'go-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.5 0],...
%     'MarkerSize',3);
% hold off
% set(findobj('color','g'),'Color',[0 0.5 0]);
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
% xlabel('Tempo (h)');
% ylabel('z (m)');
% xlim([0 120/3600])
% grid on
% 
% % Comparação entre as posições nos instantes finais
% figure(7)
% subplot(3,1,1)
% plot(tm(1:end-2)/3600,xyz_hiperdim(1:end-1,1),'bs-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.4470 0.7410],...
%     'MarkerSize',3)
% hold on
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,1),'rd-','LineWidth',1,...
%     'MarkerFaceColor','r',...
%     'MarkerSize',3)
% hold on
% plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,1),'go-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.5 0],...
%     'MarkerSize',3)
% hold off
% set(findobj('color','g'),'Color',[0 0.5 0]);
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
% xlabel('Tempo (h)');
% ylabel('x (m)');
% xlim([(tm(end)-120)/3600 tm(end)/3600])
% grid on
% legend('Medido','Estimado','Verdade', 'Location', 'southeast');
% lgd = legend;
% lgd.NumColumns = 3;
% subplot(3,1,2)
% plot(tm(1:end-2)/3600,xyz_hiperdim(1:end-1,2),'bs-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.4470 0.7410],...
%     'MarkerSize',3);
% hold on
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,2),'rd-','LineWidth',1,...
%     'MarkerFaceColor','r',...
%     'MarkerSize',3);
% hold on
% plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,2),'go-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.5 0],...
%     'MarkerSize',3);
% hold off
% set(findobj('color','g'),'Color',[0 0.5 0]);
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
% xlabel('Tempo (h)');
% ylabel('y (m)');
% xlim([(tm(end)-120)/3600 tm(end)/3600])
% grid on
% subplot(3,1,3)
% plot(tm(1:end-2)/3600,xyz_hiperdim(1:end-1,3),'bs-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.4470 0.7410],...
%     'MarkerSize',3);
% hold on
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,3),'rd-','LineWidth',1,...
%     'MarkerFaceColor','r',...
%     'MarkerSize',3);
% hold on
% plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,3),'go-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.5 0],...
%     'MarkerSize',3);
% hold off
% set(findobj('color','g'),'Color',[0 0.5 0]);
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
% xlabel('Tempo (h)');
% ylabel('z (m)');
% xlim([(tm(end)-120)/3600 tm(end)/3600])
% grid on
% 
% % Comparação entre as posições no período completo
% figure(8)
% subplot(3,1,1)
% plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,1),'g',...
%      tm(1:end-2)/3600,xyz_hiperdim(1:end-1,1),'b',...
%      tm(1:end-2)/3600,estado_aumentado(2:end,1),'r',...
%      'LineWidth',1);
% set(findobj('color','g'),'Color',[0 0.5 0]);
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
% xlabel('Tempo (h)');
% ylabel('x (m)');
% legend('Verdade','Medido','Estimado');
% grid on
% subplot(3,1,2)
% plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,2),'g',...
%      tm(1:end-2)/3600,xyz_hiperdim(1:end-1,2),'b',...
%      tm(1:end-2)/3600,estado_aumentado(2:end,2),'r',...
%     'LineWidth',1);
% set(findobj('color','g'),'Color',[0 0.5 0]);
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
% xlabel('Tempo (h)');
% ylabel('y (m)');
% grid on
% subplot(3,1,3)
% plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,3),'g',...
%      tm(1:end-2)/3600,xyz_hiperdim(1:end-1,3),'b',...
%      tm(1:end-2)/3600,estado_aumentado(2:end,3),'r',...
%     'LineWidth',1);
% set(findobj('color','g'),'Color',[0 0.5 0]);
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);xlabel('Tempo (h)');
% ylabel('z (m)');
% grid on
% 
% % Sequência de inovação
% figure(9)
% subplot(3,1,1)
% plot(tabt/3600 ,tab2(:,1),'b');
% hold on
% plot(tabt/3600,tab4(:,1),'r',tabt/3600,tab4(:,4),'r','LineWidth',1);
% hold off
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
% xlabel('Tempo (s)');
% ylabel('Sequência em x (m)');
% ylim([-80 80])
% legend('Sequência de Inovação','\pm3\sigma');
% lgd = legend;
% lgd.NumColumns = 3;
% subplot(3,1,2)
% plot(tabt/3600 ,tab2(:,2),'b');
% hold on
% plot(tabt/3600,tab4(:,1),'r',tabt/3600,tab4(:,4),'r','LineWidth',1);
% hold off
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
% xlabel('Tempo (s)');
% ylabel('Sequência em y (m)');
% ylim([-80 80])
% subplot(3,1,3)
% plot(tabt/3600 ,tab2(:,3),'b')
% hold on
% plot(tabt/3600,tab4(:,1),'r',tabt/3600,tab4(:,4),'r','LineWidth',1);
% hold off
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
% xlabel('Tempo (s)');
% ylabel('Sequência em z (m)');
% ylim([-80 80])
% 
% % Comparação de velocidades 
% figure(12)
% subplot(3,3,1)
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,4),'r',...
%     tm(1:end-2)/3600,xyz_verdade(1:end-1,4),'g',...
%      'LineWidth',1);
% title('Periodo Total')
% set(findobj('color','g'),'Color',[0 0.5 0]);
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
% xlabel('Tempo (h)');
% ylabel('v_x (m/s)');
% legend('Medido','Verdade', 'Location', 'southeast');
% lgd = legend;
% lgd.NumColumns = 3;
% grid on
% subplot(3,3,4)
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,5),'r',...
%      tm(1:end-2)/3600,xyz_verdade(1:end-1,5),'g',...
%      'LineWidth',1);
% set(findobj('color','g'),'Color',[0 0.5 0]);
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
% xlabel('Tempo (h)');
% ylabel('v_y (m/s)');
% grid on
% subplot(3,3,7)
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,6),'r',...
%      tm(1:end-2)/3600,xyz_verdade(1:end-1,6),'g',...
%      'LineWidth',1);
% set(findobj('color','g'),'Color',[0 0.5 0]);
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);xlabel('Tempo (h)');
% ylabel('v_z (m/s)');
% grid on
% subplot(3,3,2)
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,4),'rd-','LineWidth',1,...
%     'MarkerFaceColor','r',...
%     'MarkerSize',3)
% hold on
% plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,4),'go-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.5 0],...
%     'MarkerSize',3)
% hold off
% title('Instantes Iniciais')
% set(findobj('color','g'),'Color',[0 0.5 0]);
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
% xlabel('Tempo (h)');
% ylabel('v_x (m/s)');
% xlim([0 120/3600])
% grid on
% legend('Medido','Verdade', 'Location', 'southeast');
% lgd = legend;
% lgd.NumColumns = 3;
% subplot(3,3,5)
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,5),'rd-','LineWidth',1,...
%     'MarkerFaceColor','r',...
%     'MarkerSize',3);
% hold on
% plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,5),'go-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.5 0],...
%     'MarkerSize',3);
% hold off
% set(findobj('color','g'),'Color',[0 0.5 0]);
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
% xlabel('Tempo (h)');
% ylabel('v_y (m/s)');
% xlim([0 120/3600])
% grid on
% subplot(3,3,8)
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,6),'rd-','LineWidth',1,...
%     'MarkerFaceColor','r',...
%     'MarkerSize',3);
% hold on
% plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,6),'go-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.5 0],...
%     'MarkerSize',3);
% hold off
% set(findobj('color','g'),'Color',[0 0.5 0]);
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
% xlabel('Tempo (h)');
% ylabel('v_z (m/s)');
% xlim([0 120/3600])
% grid on
% subplot(3,3,3)
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,4),'rd-','LineWidth',1,...
%     'MarkerFaceColor','r',...
%     'MarkerSize',3)
% hold on
% plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,4),'go-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.5 0],...
%     'MarkerSize',3)
% title('Instantes Finais')
% hold off
% set(findobj('color','g'),'Color',[0 0.5 0]);
% set(findobj('color','b'),'Color',[0 0.4470 2.741]);
% xlabel('Tempo (h)');
% ylabel('v_x (m/s)');
% xlim([(tm(end)-120)/3600 tm(end)/3600])
% grid on
% legend('Medido','Verdade', 'Location', 'southeast');
% lgd = legend;
% lgd.NumColumns = 3;
% subplot(3,3,6)
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,5),'rd-','LineWidth',1,...
%     'MarkerFaceColor','r',...
%     'MarkerSize',3);
% hold on
% plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,5),'go-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.5 0],...
%     'MarkerSize',3);
% hold off
% set(findobj('color','g'),'Color',[0 0.5 0]);
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
% xlabel('Tempo (h)');
% ylabel('v_y (m/s)');
% xlim([(tm(end)-120)/3600 tm(end)/3600])
% grid on
% subplot(3,3,9)
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,6),'rd-','LineWidth',1,...
%     'MarkerFaceColor','r',...
%     'MarkerSize',3);
% hold on
% plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,6),'go-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.5 0],...
%     'MarkerSize',3);
% hold off
% set(findobj('color','g'),'Color',[0 0.5 0]);
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
% xlabel('Tempo (h)');
% ylabel('v_z (m/s)');
% xlim([(tm(end)-120)/3600 tm(end)/3600])
% grid on
% 
% % Comparação de posições
% figure(13)
% 
% subplot(3,3,1)
% plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,1),'g',...
%      tm(1:end-2)/3600,xyz_hiperdim(1:end-1,1),'b',...
%      tm(1:end-2)/3600,estado_aumentado(2:end,1),'r',...
%      'LineWidth',1);
% title('Periodo Total')
% set(findobj('color','g'),'Color',[0 0.5 0]);
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
% xlabel('Tempo (h)');
% ylabel('x (m)');
% legend('Verdade','Medido','Estimado', 'Location', 'southeast');
% lgd = legend;
% lgd.NumColumns = 3;
% grid on
% subplot(3,3,4)
% plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,2),'g',...
%      tm(1:end-2)/3600,xyz_hiperdim(1:end-1,2),'b',...
%      tm(1:end-2)/3600,estado_aumentado(2:end,2),'r',...
%     'LineWidth',1);
% set(findobj('color','g'),'Color',[0 0.5 0]);
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
% xlabel('Tempo (h)');
% ylabel('y (m)');
% grid on
% subplot(3,3,7)
% plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,3),'g',...
%      tm(1:end-2)/3600,xyz_hiperdim(1:end-1,3),'b',...
%      tm(1:end-2)/3600,estado_aumentado(2:end,3),'r',...
%     'LineWidth',1);
% set(findobj('color','g'),'Color',[0 0.5 0]);
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);xlabel('Tempo (h)');
% ylabel('z (m)');
% grid on
% 
% subplot(3,3,2)
% plot(tm(1:end-2)/3600,xyz_hiperdim(1:end-1,1),'bs-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.4470 0.7410],...
%     'MarkerSize',3)
% hold on
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,1),'rd-','LineWidth',1,...
%     'MarkerFaceColor','r',...
%     'MarkerSize',3)
% hold on
% plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,1),'go-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.5 0],...
%     'MarkerSize',3)
% hold off
% set(findobj('color','g'),'Color',[0 0.5 0]);
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
% xlabel('Tempo (h)');
% ylabel('x (m)');
% xlim([0 120/3600])
% title('Instantes Iniciais')
% grid on
% legend('Medido','Estimado','Verdade', 'Location', 'southeast');
% lgd = legend;
% lgd.NumColumns = 3;
% subplot(3,3,5)
% plot(tm(1:end-2)/3600,xyz_hiperdim(1:end-1,2),'bs-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.4470 0.7410],...
%     'MarkerSize',3);
% hold on
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,2),'rd-','LineWidth',1,...
%     'MarkerFaceColor','r',...
%     'MarkerSize',3);
% hold on
% plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,2),'go-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.5 0],...
%     'MarkerSize',3);
% hold off
% set(findobj('color','g'),'Color',[0 0.5 0]);
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
% xlabel('Tempo (h)');
% ylabel('y (m)');
% xlim([0 120/3600])
% grid on
% subplot(3,3,8)
% plot(tm(1:end-2)/3600,xyz_hiperdim(1:end-1,3),'bs-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.4470 0.7410],...
%     'MarkerSize',3);
% hold on
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,3),'rd-','LineWidth',1,...
%     'MarkerFaceColor','r',...
%     'MarkerSize',3);
% hold on
% plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,3),'go-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.5 0],...
%     'MarkerSize',3);
% hold off
% set(findobj('color','g'),'Color',[0 0.5 0]);
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
% xlabel('Tempo (h)');
% ylabel('z (m)');
% xlim([0 120/3600])
% grid on
% 
% subplot(3,3,3)
% plot(tm(1:end-2)/3600,xyz_hiperdim(1:end-1,1),'bs-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.4470 0.7410],...
%     'MarkerSize',3)
% hold on
% title('Instantes Finais')
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,1),'rd-','LineWidth',1,...
%     'MarkerFaceColor','r',...
%     'MarkerSize',3)
% hold on
% plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,1),'go-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.5 0],...
%     'MarkerSize',3)
% hold off
% set(findobj('color','g'),'Color',[0 0.5 0]);
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
% xlabel('Tempo (h)');
% ylabel('x (m)');
% xlim([(tm(end)-120)/3600 tm(end)/3600])
% grid on
% legend('Medido','Estimado','Verdade', 'Location', 'southeast');
% lgd = legend;
% lgd.NumColumns = 3;
% subplot(3,3,6)
% plot(tm(1:end-2)/3600,xyz_hiperdim(1:end-1,2),'bs-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.4470 0.7410],...
%     'MarkerSize',3);
% hold on
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,2),'rd-','LineWidth',1,...
%     'MarkerFaceColor','r',...
%     'MarkerSize',3);
% hold on
% plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,2),'go-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.5 0],...
%     'MarkerSize',3);
% hold off
% set(findobj('color','g'),'Color',[0 0.5 0]);
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
% xlabel('Tempo (h)');
% ylabel('y (m)');
% xlim([(tm(end)-120)/3600 tm(end)/3600])
% grid on
% subplot(3,3,9)
% plot(tm(1:end-2)/3600,xyz_hiperdim(1:end-1,3),'bs-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.4470 0.7410],...
%     'MarkerSize',3);
% hold on
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,3),'rd-','LineWidth',1,...
%     'MarkerFaceColor','r',...
%     'MarkerSize',3);
% hold on
% plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,3),'go-','LineWidth',1,...
%     'MarkerFaceColor',[0 0.5 0],...
%     'MarkerSize',3);
% hold off
% set(findobj('color','g'),'Color',[0 0.5 0]);
% set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
% xlabel('Tempo (h)');
% ylabel('z (m)');
% xlim([(tm(end)-120)/3600 tm(end)/3600])
% grid on
% 
% % Comparação entre as posições e variações de posições do lider 
% figure(14)
% subplot(2,3,1)
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,7),'r',...
%      tm(1:end-2)/3600,xyz_verdade(1:end-1,7),'g',...
%      'LineWidth',1);
% set(findobj('color','g'),'Color',[0 0.5 0]);
% xlabel('Tempo (h)');
% ylabel('$r_0$(m)', 'Interpreter', 'latex');
% legend('Estimado','Verdade');
% title('Periodo Total')
% grid on
% subplot(2,3,4)
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,8),'r',...
%     tm(1:end-2)/3600,xyz_verdade(1:end-1,8),'g',...
%     'LineWidth',1);
% set(findobj('color','g'),'Color',[0 0.5 0]);
% xlabel('Tempo (h)');
% ylabel('$\dot{r}_0$(m/s)', 'Interpreter', 'latex');
% grid on
% subplot(2,3,2)
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,7),'r',...
%      tm(1:end-2)/3600,xyz_verdade(1:end-1,7),'g',...
%      'LineWidth',1);
% set(findobj('color','g'),'Color',[0 0.5 0]);
% xlabel('Tempo (h)');
% ylabel('$r_0$(m)', 'Interpreter', 'latex');
% legend('Estimado','Verdade');
% xlim([0 120/3600])
% title('Instantes Iniciais')
% grid on
% subplot(2,3,5)
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,8),'r',...
%     tm(1:end-2)/3600,xyz_verdade(1:end-1,8),'g',...
%     'LineWidth',1);
% set(findobj('color','g'),'Color',[0 0.5 0]);
% xlabel('Tempo (h)');
% ylabel('$\dot{r}_0$(m/s)', 'Interpreter', 'latex');
% xlim([0 120/3600])
% grid on
% subplot(2,3,3)
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,7),'r',...
%      tm(1:end-2)/3600,xyz_verdade(1:end-1,7),'g',...
%      'LineWidth',1);
% set(findobj('color','g'),'Color',[0 0.5 0]);
% xlabel('Tempo (h)');
% ylabel('$r_0$(m)', 'Interpreter', 'latex');
% legend('Estimado','Verdade');
% xlim([(tm(end)-120)/3600 tm(end)/3600])
% title('Instantes Finais')
% grid on
% subplot(2,3,6)
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,8),'r',...
%     tm(1:end-2)/3600,xyz_verdade(1:end-1,8),'g',...
%     'LineWidth',1);
% set(findobj('color','g'),'Color',[0 0.5 0]);
% xlabel('Tempo (h)');
% ylabel('$\dot{r}_0$(m/s)', 'Interpreter', 'latex');
% xlim([(tm(end)-120)/3600 tm(end)/3600])
% grid on
% 
% % Comparação entre os argumentos de latitude e variações do argumento
% figure(15)
% subplot(2,3,1)
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,9),'r',...
%     tm(1:end-2)/3600,xyz_verdade(1:end-1,9),'g',...
%      'LineWidth',1);
% set(findobj('color','g'),'Color',[0 0.5 0]);
% xlabel('Tempo (h)');
% ylabel('$\theta_0$ (m)', 'Interpreter', 'latex');
% legend('Estimado','Verdade');
% title('Periodo Total')
% grid on
% subplot(2,3,4)
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,10),'r',...
%     tm(1:end-2)/3600,xyz_verdade(1:end-1,10),'g',...
%     'LineWidth',1);
% set(findobj('color','g'),'Color',[0 0.5 0]);
% xlabel('Tempo (h)');
% ylabel('$\dot{\theta}_0$ (m/s)', 'Interpreter', 'latex');
% grid on
% subplot(2,3,2)
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,9),'r',...
%     tm(1:end-2)/3600,xyz_verdade(1:end-1,9),'g',...
%      'LineWidth',1);
% set(findobj('color','g'),'Color',[0 0.5 0]);
% xlabel('Tempo (h)');
% ylabel('$\theta_0$ (m)', 'Interpreter', 'latex');
% legend('Estimado','Verdade');
% xlim([0 120/3600])
% title('Instantes Iniciais')
% grid on
% subplot(2,3,5)
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,10),'r',...
%     tm(1:end-2)/3600,xyz_verdade(1:end-1,10),'g',...
%     'LineWidth',1);
% set(findobj('color','g'),'Color',[0 0.5 0]);
% xlabel('Tempo (h)');
% ylabel('$\dot{\theta}_0$ (m/s)', 'Interpreter', 'latex');
% xlim([0 120/3600])
% grid on
% subplot(2,3,3)
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,9),'r',...
%     tm(1:end-2)/3600,xyz_verdade(1:end-1,9),'g',...
%      'LineWidth',1);
% set(findobj('color','g'),'Color',[0 0.5 0]);
% xlabel('Tempo (h)');
% ylabel('$\theta_0$ (m)', 'Interpreter', 'latex');
% legend('Estimado','Verdade');
% xlim([(tm(end)-120)/3600 tm(end)/3600])
% title('Instantes Finais')
% grid on
% subplot(2,3,6)
% plot(tm(1:end-2)/3600,estado_aumentado(2:end,10),'r',...
%     tm(1:end-2)/3600,xyz_verdade(1:end-1,10),'g',...
%     'LineWidth',1);
% set(findobj('color','g'),'Color',[0 0.5 0]);
% xlabel('Tempo (h)');
% ylabel('$\dot{\theta}_0$ (m/s)', 'Interpreter', 'latex');
% xlim([(tm(end)-120)/3600 tm(end)/3600])
% grid on
% 


%% ARQUIVO FINAL COM OS RESULTADOS
residuoxmedido = tab2(500:end,1);
residuoymedido = tab2(500:end,2);
residuozmedido = tab2(500:end,3);
residuoxestimado = tab2(500:end,4);
residuoyestimado = tab2(500:end,5);
residuozestimado = tab2(500:end,6);

rmsxmedido = rms(residuoxmedido,1); 
rmsymedido = rms(residuoymedido,1);
rmszmedido = rms(residuozmedido,1); 
rmsxestimado = rms(residuoxestimado,1);
rmsyestimado = rms(residuoyestimado,1);
rmszestimado = rms(residuozestimado,1);
por_media_x = 100*rmsxestimado/rmsxmedido;
por_media_y = 100*rmsyestimado/rmsymedido;
por_media_z = 100*rmszestimado/rmszmedido;

stdxmedido = std(residuoxmedido,1); 
stdymedido = std(residuoymedido,1);
stdzmedido = std(residuozmedido,1); 
stdxestimado = std(residuoxestimado,1);
stdyestimado = std(residuoyestimado,1);
stdzestimado = std(residuozestimado,1);
por_des_x = 100*stdxestimado/stdxmedido;
por_des_y = 100*stdyestimado/stdymedido;
por_des_z = 100*stdzestimado/stdzmedido;

diary('RESULTADOS_FINAIS.txt');
fprintf('RMSE posição medida\nx: %.10f m\ny: %.10f m\nz: %.10f m\n\n', rmsxmedido,rmsymedido,rmszmedido);
fprintf('RMSE posição estimada\nx: %.10f m\ny: %.10f m\nz: %.10f m\n\n', rmsxestimado,rmsyestimado,rmszestimado);
fprintf('Porcentagem de redução \nx: %.10f m\ny: %.10f m\nz: %.10f m\n\n\n\n',100-por_media_x,100-por_media_y,100-por_media_z);
fprintf('STD posição medida\nx: %.10f m\ny: %.10f m\nz: %.10f m\n\n', stdxmedido,stdymedido,stdzmedido);
fprintf('STD posição estimada\nx: %.10f m\ny: %.10f m\nz: %.10f m\n\n', stdxestimado,stdyestimado,stdzestimado);
fprintf('Porcentagem de redução \nx: %.10f m\ny: %.10f m\nz: %.10f m\n\n',100-por_des_x,100-por_des_y,100-por_des_z);
diary off;
