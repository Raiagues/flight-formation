% FILTRO DE KALMAN APLICADO A MANOBRAS AEROASSISTIDAS
%
format long
load data4.dat
xyz = data4(:,2:4);
%t = data1(:,1);
passo  = 5; % intervalo de amostragem das medidas (s)
%
h = size(xyz);
tabt = zeros(h(1),1);
tab1 = zeros(h(1),8);
tab2 = zeros(h(1),6);
tab3 = zeros(h(1),8);
%
% Condição inicial
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
R = diag([5e-4^2 5e-4^2 5e-4^2]); % erro da observação
%
Xo = [x1o; x2o; x3o; x4o; x5o; x6o; x7o; x8o];
Po = diag([p1o^2; p2o^2; p3o^2; p4o^2; p5o^2; p6o^2; p7o^2; p8o^2]);
%
xb = Xo;
Pb = Po;
%
H = [1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 0 0 1 0 0 0 0 0];
%
for i = 1:h(1),
    %
    % Fase de atualização
    K = (Pb*H')/(H*Pb*H' + R);
    Pc = (eye(8,8) - K*H)*Pb;
    y = [xyz(i,1); xyz(i,2); xyz(i,3)];
    res = y - H*xb;
    Xc = xb + K*res;
    res2 = Xc(1:3) - xb(1:3);
    dvp = sqrt(diag(Pc));
    %
    % Re-alimentação do vetor condição inicial
    CI = [Xc; Pc(:,1); Pc(:,2); Pc(:,3); Pc(:,4); Pc(:,5); Pc(:,6);...
        Pc(:,7); Pc(:,8)];
    %
    % Fase de propagação
    [t,Xp] = ode45(@propag,[(passo*i - passo) i*passo],CI);
    d = size(Xp);
    Xb = Xp(d(1),:)'; % estado completo propagado
    xb = Xb(1:8); % estado dinâmico propagado
    Pb = [Xb(9:16) Xb(17:24) Xb(25:32) Xb(33:40) Xb(41:48) Xb(49:56)...
    Xb(57:64) Xb(65:72)];
    %
    % Armazenamento de resultados
    tabt(i,:) = i*passo;
    tab1(i,:) = [Xc(1),Xc(2),Xc(3),Xc(4),Xc(5),Xc(6),Xc(7),Xc(8)];
    tab2(i,:) = [res',res2'];
    tab3(i,:) = [dvp(1),dvp(2),dvp(3),dvp(4),dvp(5),dvp(6),dvp(7),dvp(8)];
    %
end
%
% Organização das matrizes para plotagem dos gráficos
tm = cat(1,0,tabt);
tabx0 = [x1o x2o x3o x4o x5o x6o x7o x8o];
tab1 = cat(1,tabx0,tab1);
tabdv = tab3;
tabp0 = [p1o p2o p3o p4o p5o p6o p7o p8o];
tab3 = cat(1,tabp0,tab3);
smax = [abs(2*tabdv(:,1)) abs(2*tabdv(:,2)) abs(2*tabdv(:,3))];
smin = -smax;
tab4 = cat(2,smax,smin);
%
% Plotagem dos gráficos de resultados
figure
plot(tm,tab1(:,1:3),'LineWidth',1);
set(gca,'FontSize',14,'XMinorTick','on','YMinorTick','on');
xlabel('Tempo (s)','FontSize',14,'FontWeight','bold');
ylabel('Posição (km)','FontSize',14,'FontWeight','bold');
legend('x','y','z');
legend('boxoff');
%
figure
plot(tm,tab1(:,4:6),'LineWidth',1);
set(gca,'FontSize',14,'XMinorTick','on','YMinorTick','on');
xlabel('Tempo (s)','FontSize',14,'FontWeight','bold');
ylabel('Velocidade (km/s)','FontSize',14,'FontWeight','bold');
legend('vx','vy','vz');
legend('boxoff');
%
figure
plot(tm,tab1(:,7),'LineWidth',1);
set(gca,'FontSize',14,'XMinorTick','on','YMinorTick','on');
xlabel('Tempo (s)','FontSize',14,'FontWeight','bold');
ylabel('Constante Gravitacional Geocentrica \mu (km^3/s^2)',...
    'FontSize',14,'FontWeight','bold');
%
figure
plot(tm,tab1(:,8),'LineWidth',1);
set(gca,'FontSize',14,'XMinorTick','on','YMinorTick','on');
xlabel('Tempo (s)','FontSize',14,'FontWeight','bold');
ylabel('Coeficiente aerodinâmico de arrasto Cd (adm)',...
    'FontSize',14,'FontWeight','bold');
%
figure
plot(tabt,tab2(:,1),tabt,tab2(:,4),tabt,tab4(:,1),tabt,tab4(:,4),...
    'LineWidth',1);
set(gca,'FontSize',14,'XMinorTick','on','YMinorTick','on');
xlabel('Tempo (s)','FontSize',14,'FontWeight','bold');
ylabel('Residuo em X (km)','FontSize',14,'FontWeight','bold');
legend('Residuo de medição','Resíduo de estimação');
legend('boxoff');
%
figure
plot(tabt,tab2(:,2),tabt,tab2(:,5),tabt,tab4(:,2),tabt,tab4(:,5),...
    'LineWidth',1);
set(gca,'FontSize',14,'XMinorTick','on','YMinorTick','on');
xlabel('Tempo (s)','FontSize',14,'FontWeight','bold');
ylabel('Residuo em Y (km)','FontSize',14,'FontWeight','bold');
legend('Residuo de medição','Resíduo de estimação');
legend('boxoff');
%
figure
plot(tabt,tab2(:,3),tabt,tab2(:,6),tabt,tab4(:,3),tabt,tab4(:,6),...
    'LineWidth',1);
set(gca,'FontSize',14,'XMinorTick','on','YMinorTick','on');
xlabel('Tempo (s)','FontSize',14,'FontWeight','bold');
ylabel('Residuo em Z (km)','FontSize',14,'FontWeight','bold');
legend('Residuo de medição','Resíduo de estimação');
legend('boxoff');
%
figure
plot(tm,tab3(:,1:3),'LineWidth',1);
set(gca,'FontSize',14,'XMinorTick','on','YMinorTick','on');
xlabel('Tempo (s)','FontSize',14,'FontWeight','bold');
ylabel('Desvio Padrão na Posição (km)','FontSize',14,'FontWeight',...
    'bold');
legend('\sigma_x','\sigma_y','\sigma_z');
legend('boxoff');
%
figure
plot(tm,tab3(:,4:6),'LineWidth',1);
set(gca,'FontSize',14,'XMinorTick','on','YMinorTick','on');
xlabel('Tempo (s)','FontSize',14,'FontWeight','bold');
ylabel('Desvio Padrão na Velocidade (km/s)','FontSize',14,...
    'FontWeight','bold');
legend('\sigma_v_x','\sigma_v_y','\sigma_v_z');
legend('boxoff');
%
figure
plot(tm,tab3(:,7),'LineWidth',1);
set(gca,'FontSize',14,'XMinorTick','on','YMinorTick','on');
xlabel('Tempo (s)','FontSize',14,'FontWeight','bold');
ylabel('Desvio Padrão na Constante Grav. Geoc. (km^3/s^2)',...
    'FontSize',14,'FontWeight','bold');
%
figure
plot(tm,tab3(:,8),'LineWidth',1);
set(gca,'FontSize',14,'XMinorTick','on','YMinorTick','on');
xlabel('Tempo (s)','FontSize',14,'FontWeight','bold');
ylabel('Desvio Padrão no Coeficiente de arrasto (adm)',...
    'FontSize',14,'FontWeight','bold');
%
figure
plot(xyz(:,1),xyz(:,2),tab1(:,1),tab1(:,2),'LineWidth',1);
set(gca,'FontSize',14,'XMinorTick','on','YMinorTick','on');
xlabel('Posição em x (km)','FontSize',14,'FontWeight','bold');
ylabel('Posição em y (km)','FontSize',14,'FontWeight','bold');
legend('Órbita medida','Órbita estimada');
legend('boxoff');
%
% Melhor ajuste de parâmetros encontrado
% x1o = 6.498139999999999e3; % (X) km
% x2o = 0; % (Y) km
% x3o = 0; % (Z) km
% x4o = 0; % (Vx) km/s
% x5o = 8.076552708163819; % (Vy) km/s
% x6o = 0; % (Vz) km/s
% x7o = 398600.4418; % (mi) km^3/s^2
% x8o = 2; % (Cd) adme3; % (X) km
% %
% p1o = 1; % (Sigma X) km
% p2o = 1; % (Sigma Y) km
% p3o = 1; % (Sigma Z) km
% p4o = 0.001; % (Sigma Vx) km/s
% p5o = 0.001; % (Sigma Vy) km/s
% p6o = 0.001; % (Sigma Vz) km/s
% p7o = 0.1; % (Sigma mi) km^3/s^2
% p8o = 0.2; % (Sigma Cd) adm
% %
% w1 = 1*1e-5;
% w2 = 1*1e-5;
% w3 = 1*1e-5;
% w4 = 1*1e-10;
% w5 = 1*1e-10;
% %
% r = 0.01;