%clc
%clear all
%% Inicialização Código
%Parâmetros gerais
global mu
mu = 3.986004418e14;
Earth_r = 6378140;
a0 = Earth_r + 600e3; 
n0 = sqrt(mu/a0^3);
e0 = 0.1;
omega0 = 0;
f0 = 0;
theta = omega0 + f0;

%Importando medidas iniciais
format long
load VALORESINICIAIS.dat
xyz = VALORESINICIAIS(:,2:4);
passo  = 5; 
h = size(xyz);

%Definindo matrizes para guardar resultados
tabt = zeros(h(1),1);
tab1 = zeros(h(1),10);
tab2 = zeros(h(1),6);
tab3 = zeros(h(1),10);

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

%Valores iniciais para vetor de estados
x1o = xn*a0;            % x [m]
x2o = yn*a0;            % y [m]
x3o = zn*a0;            % z [m]
x4o = xdot_n*(n0*a0);   % xdot [m/s]
x5o = ydot_n*(n0*a0);   % ydot [m/s]
x6o = zdot_n*(n0*a0);   % zdot [m/s]
x7o = r0;               % r0 [m]
x8o = r0dot;            % r0dot [m]
x9o = theta;            % theta0 [rad]
x10o = thetadot;        % theta0dot [rad/s]

%Valores iniciais para matriz de covariancia
p1o = 1000;              % sigma_x [m]
p2o = 1000;              % sigma_y [m]
p3o = 1000;              % sigma_z [m]
p4o = 100;              % sigma_xdot [m/s]
p5o = 100;              % sigma_ydot [m/s]
p6o = 100;              % sigma_zdot [m/s]
p7o = 0.1;                % sigma_r0 [m]
p8o = 0.1;                % sigma_r0dot [m]
p9o = 0.1;                % sigma_theta0 [rad]
p10o = 0.1;               % sigma_theta0dot [rad/s]

%% Inicialização do Filtro de Kalman
%Vetores de estado e ruído
Xo = [x1o; x2o; x3o; x4o; x5o; x6o; x7o; x8o; x9o; x10o];
Po = diag([p1o^2; p2o^2; p3o^2; p4o^2; p5o^2; p6o^2; p7o^2; p8o^2; p9o^2; p10o^2]);

%Matriz erro de observação
R = diag([5e-4^2 5e-4^2 5e-4^2]);

%Matriz condições iniciais
Xc = Xo;
Pc = Po;
%CI = [Xc; Pc(:,1); Pc(:,2); Pc(:,3); Pc(:,4); Pc(:,5); Pc(:,6); Pc(:,7); Pc(:,8); Pc(:,9); Pc(:,10)];

%Matriz de observação
H = [1 0 0 0 0 0 0 0 0 0;0 1 0 0 0 0 0 0 0 0;0 0 1 0 0 0 0 0 0 0];

%% Filtro de Kalman
for i = 1:h(1),
    %Fase Atualização
    K = Pc*H'*inv(H*Pc*H'+R);
    Pc = (eye(10,10)-K*H)*Pc;
    y = [xyz(i,1); xyz(i,2); xyz(i,3)];
    
    %Calculando resíduos
    res = y-H*xb;
    Xc = xb+K*res;
    res2 = Xc(1:3)-xb(1:3);
    dvp = sqrt(diag(Pc));
    
    %Atualizando condições iniciais
    CI = [Xc; Pc(:,1); Pc(:,2); Pc(:,3); Pc(:,4); Pc(:,5); Pc(:,6); Pc(:,7); Pc(:,8); Pc(:,9); Pc(:,10)];
    
    %Fase Propagação
    [t,Xp] = ode45(@edo,[(passo*i-passo) i*passo],CI);
    d = size(Xp);
    Xb = Xp(d(1),:)';
    xb = Xb(1:10);
    Pb = [Xb(11:20) Xb(21:30) Xb(31:40) Xb(41:50) Xb(51:60) Xb(61:70) Xb(71:80) Xb(81:90) Xb(91:100) Xb(101:110)];
    

    %Guardando variáveis
    tabt(i,:) = i*passo;
    tab1(i,:) = [Xc(1),Xc(2),Xc(3),Xc(4),Xc(5),Xc(6),Xc(7),Xc(8),Xc(9),Xc(10)];
    tab2(i,:) = [res',res2'];
    tab3(i,:) = [dvp(1),dvp(2),dvp(3),dvp(4),dvp(5),dvp(6),dvp(7),dvp(8),dvp(9),dvp(10)];
end

%% Organizando dados para plotagem
tm = cat(1,0,tabt);
tabx0 = [x1o x2o x3o x4o x5o x6o x7o x8o x9o x10o];
tab1 = cat(1,tabx0,tab1);
tabdv = tab3;
tabp0 = [p1o p2o p3o p4o p5o p6o p7o p8o p9o p10o];
tab3 = cat(1,tabp0,tab3);
smax = [abs(2*tabdv(:,1)) abs(2*tabdv(:,2)) abs(2*tabdv(:,3))];
smin = -smax;
tab4 = cat(2,smax,smin);

%% Plotagens
%Comparação entre órbita medida e órbita estimada
figure(1)
subplot(2,2,1)
plot(xyz(:,1)/a0,xyz(:,2)/a0)
hold on
plot(tab1(:,1)/a0,tab1(:,2)/a0)
xlim([-0.06 0.06]);
ylim([-0.08 0.08]);
xlabel('x/a0');
ylabel('y/a0');
legend('Órbita Medida','Órbita Estimada','Location','NorthWest');
grid on

subplot(2,2,2)
plot(xyz(:,1)/a0,xyz(:,3)/a0)
hold on
plot(tab1(:,1)/a0,tab1(:,3)/a0)
xlim([-0.12 0.12]);
ylim([-0.14 0.14]);
xlabel('x/a0');
ylabel('z/a0');
grid on

subplot(2,2,3)
plot(xyz(:,2)/a0,xyz(:,3)/a0)
hold on
plot(tab1(:,2)/a0,tab1(:,3)/a0)
xlim([-0.12 0.12]);
ylim([-0.14 0.14]);
xlabel('y/a0');
ylabel('z/a0');
grid on

subplot(2,2,4)
plot3(xyz(:,1)/a0,xyz(:,2)/a0,xyz(:,3)/a0)
hold on
plot3(tab1(:,1)/a0,tab1(:,2)/a0,tab1(:,3)/a0)
xlim([-0.12 0.12]);
ylim([-0.12 0.12]);
zlim([-0.12 0.12]);
xlabel('x/a0');
ylabel('y/a0');
zlabel('z/a0');
hold on
plot3(xyz(1,1)/a0,xyz(1,2)/a0,xyz(1,3)/a0,'*r')
grid on


figure(6)
subplot(3,1,1)
plot(tabt,tab2(:,1),'b',tabt,tab2(:,4),'r',tabt,tab4(:,1),'y',tabt,tab4(:,4),'y','LineWidth',1);
set(findobj('color','g'),'Color',[0 0.5 0]);
set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
set(findobj('color','y'),'Color',[0.9290 0.6940 0.1250]);
xlabel('Tempo (s)');
ylabel('Residuo em x (m)');
xlim([0 4000])
legend('Residuo de medição','Resíduo de estimação','\pm\sigma');
lgd = legend;
lgd.NumColumns = 3;
subplot(3,1,2)
plot(tabt,tab2(:,2),'b',tabt,tab2(:,5),'r',tabt,tab4(:,2),'y',tabt,tab4(:,5),'y','LineWidth',1);
set(findobj('color','g'),'Color',[0 0.5 0]);
set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
set(findobj('color','y'),'Color',[0.9290 0.6940 0.1250]);
xlabel('Tempo (s)');
ylabel('Residuo em y (m)');
xlim([0 4000])
subplot(3,1,3)
plot(tabt,tab2(:,3),'b',tabt,tab2(:,6),'r',tabt,tab4(:,3),'y',tabt,tab4(:,6),'y','LineWidth',1);
set(findobj('color','g'),'Color',[0 0.5 0]);
set(findobj('color','b'),'Color',[0 0.4470 0.7410]);
set(findobj('color','y'),'Color',[0.9290 0.6940 0.1250]);
xlabel('Tempo (s)');
ylabel('Residuo em z (m)');
xlim([0 4000])
