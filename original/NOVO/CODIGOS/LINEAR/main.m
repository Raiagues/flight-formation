clc
clear all
format long
%% Inicialização Código
%Parâmetros gerais
global mu                           % Definindo variável global
mu = 3.986004418e14;                % Parametro gravitacional da terra
Earth_r = 6378140;                  % Raio da Terra
a0 = Earth_r + 600e3;               % Distância do Líder
n0 = sqrt(mu/a0^3);                 %
T = 2*pi*sqrt(a0^3/mu);             % Período

% Importando medidas iniciais
format long
load data.dat               % Coordenadas medidas
load data_verdade.dat       % Coordenadas sem ruídos
load dados_hiperdim.dat     % Coordenadas com ruídos hiperdimensionados

% Declarando coordenadas sem ruídos
xyz_verdade = data_verdade(:,2:7);         % Coordenadas sem ruídos
xyz = data(:,2:4);                          % Coordenadas medidas
xyz_hiperdim = dados_hiperdim(:,2:4);       % Coordenadas com ruídos hiperdimensionados
passo  = 5;                                 % Passo de amostragem
h = size(xyz);

tabt = zeros(h(1),1);
tab1 = zeros(h(1),6);
tab2 = zeros(h(1),6);
tab3 = zeros(h(1),6);
  
%% Condições iniciais
x1o = 69.78e3;          % x                 [m]
x2o = 139.56e3;         % y                 [m]
x3o = 104.67e3;         % z                 [m]
x4o = 7.5579;           % xdot              [m/s]
x5o = -151.116;         % ydot              [m/s]
x6o = 15.116;           % zdot              [m/s]

% Valores iniciais para matriz de covariancia
p1o = 20;               % sigma_x            [m]
p2o = 20;               % sigma_y            [m]
p3o = 20;               % sigma_z            [m]
p4o = 0.1;              % sigma_vx           [m/s]
p5o = 0.1;              % sigma_vy           [m/s]
p6o = 0.1;              % sigma_vz           [m/s]

%% Inicialização do Filtro de Kalman
% Vetor de estado e matriz de covariância
Xo = [x1o; x2o; x3o; x4o; x5o; x6o];
Po = diag([p1o^2; p2o^2; p3o^2; p4o^2; p5o^2; p6o^2]);

% Matriz erro de observação
R = diag([(10*sqrt(2))^2 (10*sqrt(2))^2 (10*sqrt(2))^2]); % [m]

% Matrizes condições iniciais
xb = Xo;                                                  % [m]
Pb = Po;                                                  % [m]

% Matriz de observação
H = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0];

for i = 1:h(1),
    % Fase Atualização
    K = Pb*H'*inv(H*Pb*H'+R);                       
    Pc = (eye(6,6)-K*H)*Pb;
    y = [xyz(i,1); xyz(i,2); xyz(i,3)];
    y_verdade = [xyz_verdade(i,1); xyz_verdade(i,2); xyz_verdade(i,3)];

    % Calculando sequência de inovação
    res = y-H*xb;
    Xc = xb+K*res;
    res2 = Xc(1:3)-xb(1:3);
    dvp = sqrt(diag(Pc));
      
    %Atualizando condições iniciais
    CI = [Xc; Pc(:,1); Pc(:,2); Pc(:,3); Pc(:,4); Pc(:,5); Pc(:,6)];
 
    % Fase Propagação 
    [t,Xp] = ode45(@edo,[(passo*i-passo) i*passo],CI);
    d = size(Xp);
    Xb = Xp(d(1),:)';
    xb = Xb(1:6);
    Pc = [Xb(7:12) Xb(13:18) Xb(19:24) Xb(25:30) Xb(31:36) Xb(37:42)];
 
   %Guardando variáveis
    tabt(i,:) = i*passo;
    tab1(i,:) = [Xc(1),Xc(2),Xc(3),Xc(4),Xc(5),Xc(6)];
    tab2(i,:) = [res',res2'];
    tab3(i,:) = [dvp(1),dvp(2),dvp(3),dvp(4),dvp(5),dvp(6)];
end

%% Organizando dados para plotagem
tm = cat(1,0,tabt);
tabx0 = [x1o x2o x3o x4o x5o x6o];
tab1 = cat(1,tabx0,tab1);
tabdv = tab3;
tabp0 = [p1o p2o p3o p4o p5o p6o];
tab3 = cat(1,tabp0,tab3);
smax = [abs(3*tabdv(:,1)) abs(3*tabdv(:,2)) abs(3*tabdv(:,3))];
smin = -smax;
tab4 = cat(2,smax,smin);

% Retirando os ruídos do vetor de medidas
resx = (tab1(2:end,1)-xyz_verdade(:,1))*1;
resy = (tab1(2:end,2)-xyz_verdade(:,2))*1;
resz = (tab1(2:end,3)-xyz_verdade(:,3))*1;
resvx = (tab1(2:end,4)-xyz_verdade(:,4));
resvy = (tab1(2:end,5)-xyz_verdade(:,5));
resvz = (tab1(2:end,6)-xyz_verdade(:,6));

% Criando um novo vetor com os ruídos hiperdimensinados
estado_aumentado = [xyz_verdade(:,1)+resx,xyz_verdade(:,2)+resy,xyz_verdade(:,3)+resz,...
                    xyz_verdade(:,4)+resvx,xyz_verdade(:,5)+resvy,xyz_verdade(:,6)+resvz];
  
%% Plotagens
% Comparação entre órbita medida e órbita estimada
figure(1)
subplot(2,2,1)
plot(xyz_hiperdim(2:end,1)/a0,xyz_hiperdim(2:end,2)/a0)
hold on
plot(estado_aumentado(2:end,1)/a0,estado_aumentado(2:end,2)/a0)
xlim([-0.02 0.02]);
ylim([-0.01 0.04]);
xlabel('x/a0');
ylabel('y/a0');
legend('Órbita Medida','Órbita Estimada','Location','NorthWest');
grid on
subplot(2,2,2)
plot(xyz_hiperdim(2:end,1)/a0,xyz_hiperdim(2:end,3)/a0)
hold on
plot(estado_aumentado(2:end,1)/a0,estado_aumentado(2:end,3)/a0)
xlim([-0.02 0.02]);
ylim([-0.02 0.02]);
xlabel('x/a0');
ylabel('z/a0');
grid on
subplot(2,2,3)
plot(xyz_hiperdim(2:end,2)/a0,xyz_hiperdim(2:end,3)/a0)
hold on
plot(estado_aumentado(2:end,2)/a0,estado_aumentado(2:end,3)/a0)
xlim([-0.02 0.04]);
ylim([-0.02 0.02]);
xlabel('y/a0');
ylabel('z/a0');
grid on
subplot(2,2,4)
plot3(xyz_hiperdim(2:end,1)/a0,xyz_hiperdim(2:end,2)/a0,xyz_hiperdim(2:end,3)/a0)
hold on
plot3(estado_aumentado(2:end,1)/a0,estado_aumentado(2:end,2)/a0,estado_aumentado(2:end,3)/a0)
xlim([-0.02 0.04]);
ylim([-0.02 0.04]);
zlim([-0.02 0.04]);
xlabel('x/a0');
ylabel('y/a0');
zlabel('z/a0');
hold on
plot3(xyz(1,1)/a0,xyz(1,2)/a0,xyz(1,3)/a0,'*r')
grid on

% Desvio padrão da posição
figure(2)
desvio_x = tab3(:,1);
desvio_y = tab3(:,2);
desvio_z = tab3(:,3);
subplot(3,1,1)
plot(tm/3600, desvio_x, 'LineWidth',1);
xlim([0 0.5])
grid on
xlabel('Tempo (h)');
ylabel('$\sigma_{x}$(m)', 'Interpreter', 'latex');
subplot(3,1,2)
plot(tm/3600, desvio_y, 'LineWidth',1);
xlim([0 0.5])
grid on
xlabel('Tempo (h)');
ylabel('$\sigma_{y}$(m)', 'Interpreter', 'latex');
subplot(3,1,3)
plot(tm/3600, desvio_z, 'LineWidth',1);
xlim([0 0.5])
grid on
xlabel('Tempo (h)');
ylabel('$\sigma_{z}$(m)', 'Interpreter', 'latex');

% Desvio padrão da velocidade
figure(3)
desvio_vx = tab3(:,4);
desvio_vy = tab3(:,5);
desvio_vz = tab3(:,6);
subplot(3,1,1)
plot(tm/3600, desvio_vx, 'LineWidth',1);
xlim([0 0.5])
grid on
xlabel('Tempo (h)');
ylabel('$\sigma_{v_x}$(m/s)', 'Interpreter', 'latex');
subplot(3,1,2)
plot(tm/3600, desvio_vy, 'LineWidth',1);
xlim([0 0.5])
grid on
xlabel('Tempo (h)');
ylabel('$\sigma_{v_y}$(m/s)', 'Interpreter', 'latex');
subplot(3,1,3)
plot(tm/3600, desvio_vz, 'LineWidth',1);
xlim([0 0.5])
grid on
xlabel('Tempo (h)');
ylabel('$\sigma_{v_0}$(m/s)', 'Interpreter', 'latex');

% Sequência de inovação
figure(4)
subplot(3,1,1)
plot(tabt/3600 ,tab2(:,1),'b');
hold on
plot(tabt/3600,tab4(:,1),'r',tabt/3600,tab4(:,4),'r','LineWidth',1);
hold off
set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
xlabel('Tempo (s)');
ylabel('Sequência em x (m)');
ylim([-80 80])
legend('Sequência de Inovação','\pm3\sigma');
lgd = legend;
lgd.NumColumns = 3;
subplot(3,1,2)
plot(tabt/3600 ,tab2(:,2),'b');
hold on
plot(tabt/3600,tab4(:,1),'r',tabt/3600,tab4(:,4),'r','LineWidth',1);
hold off
set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
xlabel('Tempo (s)');
ylabel('Sequência em y (m)');
ylim([-80 80])
subplot(3,1,3)
plot(tabt/3600 ,tab2(:,3),'b')
hold on
plot(tabt/3600,tab4(:,1),'r',tabt/3600,tab4(:,4),'r','LineWidth',1);
hold off
set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
xlabel('Tempo (s)');
ylabel('Sequência em z (m)');
ylim([-80 80])

% Comparação de velocidades 
figure(5)
subplot(3,3,1)
plot(tm(1:end-2)/3600,estado_aumentado(2:end,4),'r',...
    tm(1:end-2)/3600,xyz_verdade(1:end-1,4),'g',...
     'LineWidth',1);
title('Periodo Total')
set(findobj('color','g'),'Color',[0 0.5 0]);
set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
xlabel('Tempo (h)');
ylabel('v_x (m/s)');
legend('Medido','Verdade', 'Location', 'southeast');
lgd = legend;
lgd.NumColumns = 3;
grid on
subplot(3,3,4)
plot(tm(1:end-2)/3600,estado_aumentado(2:end,5),'r',...
     tm(1:end-2)/3600,xyz_verdade(1:end-1,5),'g',...
     'LineWidth',1);
set(findobj('color','g'),'Color',[0 0.5 0]);
set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
xlabel('Tempo (h)');
ylabel('v_y (m/s)');
grid on
subplot(3,3,7)
plot(tm(1:end-2)/3600,estado_aumentado(2:end,6),'r',...
     tm(1:end-2)/3600,xyz_verdade(1:end-1,6),'g',...
     'LineWidth',1);
set(findobj('color','g'),'Color',[0 0.5 0]);
set(findobj('color','b'),'Color',[0 0.4470 0.7410]);xlabel('Tempo (h)');
ylabel('v_z (m/s)');
grid on
subplot(3,3,2)
plot(tm(1:end-2)/3600,estado_aumentado(2:end,4),'rd-','LineWidth',1,...
    'MarkerFaceColor','r',...
    'MarkerSize',3)
hold on
plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,4),'go-','LineWidth',1,...
    'MarkerFaceColor',[0 0.5 0],...
    'MarkerSize',3)
hold off
title('Instantes Iniciais')
set(findobj('color','g'),'Color',[0 0.5 0]);
set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
xlabel('Tempo (h)');
ylabel('v_x (m/s)');
xlim([0 120/3600])
grid on
legend('Medido','Verdade', 'Location', 'southeast');
lgd = legend;
lgd.NumColumns = 3;
subplot(3,3,5)
plot(tm(1:end-2)/3600,estado_aumentado(2:end,5),'rd-','LineWidth',1,...
    'MarkerFaceColor','r',...
    'MarkerSize',3);
hold on
plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,5),'go-','LineWidth',1,...
    'MarkerFaceColor',[0 0.5 0],...
    'MarkerSize',3);
hold off
set(findobj('color','g'),'Color',[0 0.5 0]);
set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
xlabel('Tempo (h)');
ylabel('v_y (m/s)');
xlim([0 120/3600])
grid on
subplot(3,3,8)
plot(tm(1:end-2)/3600,estado_aumentado(2:end,6),'rd-','LineWidth',1,...
    'MarkerFaceColor','r',...
    'MarkerSize',3);
hold on
plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,6),'go-','LineWidth',1,...
    'MarkerFaceColor',[0 0.5 0],...
    'MarkerSize',3);
hold off
set(findobj('color','g'),'Color',[0 0.5 0]);
set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
xlabel('Tempo (h)');
ylabel('v_z (m/s)');
xlim([0 120/3600])
grid on
subplot(3,3,3)
plot(tm(1:end-2)/3600,estado_aumentado(2:end,4),'rd-','LineWidth',1,...
    'MarkerFaceColor','r',...
    'MarkerSize',3)
hold on
plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,4),'go-','LineWidth',1,...
    'MarkerFaceColor',[0 0.5 0],...
    'MarkerSize',3)
title('Instantes Finais')
hold off
set(findobj('color','g'),'Color',[0 0.5 0]);
set(findobj('color','b'),'Color',[0 0.4470 2.741]);
xlabel('Tempo (h)');
ylabel('v_x (m/s)');
xlim([(tm(end)-120)/3600 tm(end)/3600])
grid on
legend('Medido','Verdade', 'Location', 'southeast');
lgd = legend;
lgd.NumColumns = 3;
subplot(3,3,6)
plot(tm(1:end-2)/3600,estado_aumentado(2:end,5),'rd-','LineWidth',1,...
    'MarkerFaceColor','r',...
    'MarkerSize',3);
hold on
plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,5),'go-','LineWidth',1,...
    'MarkerFaceColor',[0 0.5 0],...
    'MarkerSize',3);
hold off
set(findobj('color','g'),'Color',[0 0.5 0]);
set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
xlabel('Tempo (h)');
ylabel('v_y (m/s)');
xlim([(tm(end)-120)/3600 tm(end)/3600])
grid on
subplot(3,3,9)
plot(tm(1:end-2)/3600,estado_aumentado(2:end,6),'rd-','LineWidth',1,...
    'MarkerFaceColor','r',...
    'MarkerSize',3);
hold on
plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,6),'go-','LineWidth',1,...
    'MarkerFaceColor',[0 0.5 0],...
    'MarkerSize',3);
hold off
set(findobj('color','g'),'Color',[0 0.5 0]);
set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
xlabel('Tempo (h)');
ylabel('v_z (m/s)');
xlim([(tm(end)-120)/3600 tm(end)/3600])
grid on

% Comparação de posições
figure(6)
subplot(3,3,1)
plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,1),'g',...
     tm(1:end-2)/3600,xyz_hiperdim(1:end-1,1),'b',...
     tm(1:end-2)/3600,estado_aumentado(1:end-1,1),'r',...
     'LineWidth',1);
title('Periodo Total')
set(findobj('color','g'),'Color',[0 0.5 0]);
set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
xlabel('Tempo (h)');
ylabel('x (m)');
legend('Verdade','Medido','Estimado', 'Location', 'southeast');
lgd = legend;
lgd.NumColumns = 3;
grid on
subplot(3,3,4)
plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,2),'g',...
     tm(1:end-2)/3600,xyz_hiperdim(1:end-1,2),'b',...
     tm(1:end-2)/3600,estado_aumentado(1:end-1,2),'r',...
    'LineWidth',1);
set(findobj('color','g'),'Color',[0 0.5 0]);
set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
xlabel('Tempo (h)');
ylabel('y (m)');
grid on
subplot(3,3,7)
plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,3),'g',...
     tm(1:end-2)/3600,xyz_hiperdim(1:end-1,3),'b',...
     tm(1:end-2)/3600,estado_aumentado(1:end-1,3),'r',...
    'LineWidth',1);
set(findobj('color','g'),'Color',[0 0.5 0]);
set(findobj('color','b'),'Color',[0 0.4470 0.7410]);xlabel('Tempo (h)');
ylabel('z (m)');
grid on

subplot(3,3,2)
plot(tm(1:end-2)/3600,xyz_hiperdim(1:end-1,1),'bs-','LineWidth',1,...
    'MarkerFaceColor',[0 0.4470 0.7410],...
    'MarkerSize',3)
hold on
plot(tm(1:end-2)/3600,estado_aumentado(1:end-1,1),'rd-','LineWidth',1,...
    'MarkerFaceColor','r',...
    'MarkerSize',3)
hold on
plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,1),'go-','LineWidth',1,...
    'MarkerFaceColor',[0 0.5 0],...
    'MarkerSize',3)
hold off
set(findobj('color','g'),'Color',[0 0.5 0]);
set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
xlabel('Tempo (h)');
ylabel('x (m)');
xlim([0 120/3600])
title('Instantes Iniciais')
grid on
legend('Medido','Estimado','Verdade', 'Location', 'southeast');
lgd = legend;
lgd.NumColumns = 3;
subplot(3,3,5)
plot(tm(1:end-2)/3600,xyz_hiperdim(1:end-1,2),'bs-','LineWidth',1,...
    'MarkerFaceColor',[0 0.4470 0.7410],...
    'MarkerSize',3);
hold on
plot(tm(1:end-2)/3600,estado_aumentado(1:end-1,2),'rd-','LineWidth',1,...
    'MarkerFaceColor','r',...
    'MarkerSize',3);
hold on
plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,2),'go-','LineWidth',1,...
    'MarkerFaceColor',[0 0.5 0],...
    'MarkerSize',3);
hold off
set(findobj('color','g'),'Color',[0 0.5 0]);
set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
xlabel('Tempo (h)');
ylabel('y (m)');
xlim([0 120/3600])
grid on
subplot(3,3,8)
plot(tm(1:end-2)/3600,xyz_hiperdim(1:end-1,3),'bs-','LineWidth',1,...
    'MarkerFaceColor',[0 0.4470 0.7410],...
    'MarkerSize',3);
hold on
plot(tm(1:end-2)/3600,estado_aumentado(1:end-1,3),'rd-','LineWidth',1,...
    'MarkerFaceColor','r',...
    'MarkerSize',3);
hold on
plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,3),'go-','LineWidth',1,...
    'MarkerFaceColor',[0 0.5 0],...
    'MarkerSize',3);
hold off
set(findobj('color','g'),'Color',[0 0.5 0]);
set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
xlabel('Tempo (h)');
ylabel('z (m)');
xlim([0 120/3600])
grid on

subplot(3,3,3)
plot(tm(1:end-2)/3600,xyz_hiperdim(1:end-1,1),'bs-','LineWidth',1,...
    'MarkerFaceColor',[0 0.4470 0.7410],...
    'MarkerSize',3)
hold on
title('Instantes Finais')
plot(tm(1:end-2)/3600,estado_aumentado(1:end-1,1),'rd-','LineWidth',1,...
    'MarkerFaceColor','r',...
    'MarkerSize',3)
hold on
plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,1),'go-','LineWidth',1,...
    'MarkerFaceColor',[0 0.5 0],...
    'MarkerSize',3)
hold off
set(findobj('color','g'),'Color',[0 0.5 0]);
set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
xlabel('Tempo (h)');
ylabel('x (m)');
xlim([(tm(end)-120)/3600 tm(end)/3600])
grid on
legend('Medido','Estimado','Verdade', 'Location', 'southeast');
lgd = legend;
lgd.NumColumns = 3;
subplot(3,3,6)
plot(tm(1:end-2)/3600,xyz_hiperdim(1:end-1,2),'bs-','LineWidth',1,...
    'MarkerFaceColor',[0 0.4470 0.7410],...
    'MarkerSize',3);
hold on
plot(tm(1:end-2)/3600,estado_aumentado(1:end-1,2),'rd-','LineWidth',1,...
    'MarkerFaceColor','r',...
    'MarkerSize',3);
hold on
plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,2),'go-','LineWidth',1,...
    'MarkerFaceColor',[0 0.5 0],...
    'MarkerSize',3);
hold off
set(findobj('color','g'),'Color',[0 0.5 0]);
set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
xlabel('Tempo (h)');
ylabel('y (m)');
xlim([(tm(end)-120)/3600 tm(end)/3600])
grid on
subplot(3,3,9)
plot(tm(1:end-2)/3600,xyz_hiperdim(1:end-1,3),'bs-','LineWidth',1,...
    'MarkerFaceColor',[0 0.4470 0.7410],...
    'MarkerSize',3);
hold on
plot(tm(1:end-2)/3600,estado_aumentado(1:end-1,3),'rd-','LineWidth',1,...
    'MarkerFaceColor','r',...
    'MarkerSize',3);
hold on
plot(tm(1:end-2)/3600,xyz_verdade(1:end-1,3),'go-','LineWidth',1,...
    'MarkerFaceColor',[0 0.5 0],...
    'MarkerSize',3);
hold off
set(findobj('color','g'),'Color',[0 0.5 0]);
set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
xlabel('Tempo (h)');
ylabel('z (m)');
xlim([(tm(end)-120)/3600 tm(end)/3600])
grid on


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

por_media_x = 100- 100*rmsxestimado/rmsxmedido;
por_media_y = 100 - 100*rmsyestimado/rmsymedido;
por_media_z = 100 - 100*rmszestimado/rmszmedido;

stdxmedido = std(residuoxmedido,1); 
stdymedido = std(residuoymedido,1);
stdzmedido = std(residuozmedido,1); 
stdxestimado = std(residuoxestimado,1);
stdyestimado = std(residuoyestimado,1);


stdzestimado = std(residuozestimado,1);
por_des_x = 100 - 100*stdxestimado/stdxmedido;
por_des_y = 100 - 100*stdyestimado/stdymedido;
por_des_z = 100 - 100*stdzestimado/stdzmedido;

diary('RESULTADOS_FINAIS.txt');
fprintf('RMSE posição medida\nx: %.10f m\ny: %.10f m\nz: %.10f m\n\n', rmsxmedido,rmsymedido,rmszmedido);
fprintf('RMSE posição estimada\nx: %.10f m\ny: %.10f m\nz: %.10f m\n\n', rmsxestimado,rmsyestimado,rmszestimado);
fprintf('Porcentagem de redução \nx: %.10f m\ny: %.10f m\nz: %.10f m\n\n\n\n',por_media_x,por_media_y,por_media_z);
fprintf('STD posição medida\nx: %.10f m\ny: %.10f m\nz: %.10f m\n\n', stdxmedido,stdymedido,stdzmedido);
fprintf('STD posição estimada\nx: %.10f m\ny: %.10f m\nz: %.10f m\n\n', stdxestimado,stdyestimado,stdzestimado);
fprintf('Porcentagem de redução \nx: %.10f m\ny: %.10f m\nz: %.10f m\n\n',por_des_x,por_des_y,por_des_z);
diary off;



figure(10)
subplot(3,1,1)
plot(tm(1:end-2)/3600,estado_aumentado(2:end,5),'r',...
    tm(1:end-2)/3600,xyz_verdade(1:end-1,5),'g',...
     'LineWidth',1);
title('Periodo Total')
set(findobj('color','g'),'Color',[0 0.5 0]);
set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
xlabel('Tempo (h)');
ylabel('v_x (m/s)');
legend('Medido','Verdade', 'Location', 'southeast');
lgd = legend;
lgd.NumColumns = 3;
grid on
subplot(3,1,2)
plot(tm(1:end-2)/3600,estado_aumentado(2:end,5),'r',...
     tm(1:end-2)/3600,xyz_verdade(1:end-1,5),'g',...
     'LineWidth',1);
set(findobj('color','g'),'Color',[0 0.5 0]);
set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
xlabel('Tempo (h)');
ylabel('v_y (m/s)');
title('Próximo à mudança de direção')

grid on
subplot(3,1,3)
plot(tm(1:end-2)/3600,estado_aumentado(2:end,5),'r',...
     tm(1:end-2)/3600,xyz_verdade(1:end-1,5),'g',...
     'LineWidth',1);
set(findobj('color','g'),'Color',[0 0.5 0]);
set(findobj('color','b'),'Color',[0 0.4470 0.7410]);xlabel('Tempo (h)');
ylabel('v_z (m/s)');
title('Afastado da mudança de direção')
grid on