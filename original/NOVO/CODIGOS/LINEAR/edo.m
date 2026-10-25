function Xp = edo(t,ci)

global mu
Earth_r = 6378140;
a = Earth_r + 600e3;
n = sqrt(mu/a^3);

q1o = 0.1;                  % q_x             [m]
q2o = 0.1;                  % q_y             [m]
q3o = 0.1;                  % q_z             [m]
q4o = 0.01;               % q_vx            [m/s]
q5o = 0.01;               % q_vy            [m/s]
q6o = 0.01;               % q_vz            [m/s]

Q = diag([q1o^2; q2o^2; q3o^2; q4o^2; q5o^2; q6o^2]);

x1 = ci(1);
x2 = ci(2);
x3 = ci(3);
x4 = ci(4);
x5 = ci(5);
x6 = ci(6);
p = [ci(7:12) ci(13:18) ci(19:24) ci(25:30) ci(31:36) ci(37:42)];

df4x1 = 3*n^2;
df4x5 = 2*n;
df5x4 = -2*n;
df6x3 = -n^2;
F = [0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1;
    df4x1 0 0 0 df4x5 0;
    0 0 0 df5x4 0 0;
    0 0 df6x3 0 0 0];

dx1 = x4;
dx2 = x5;
dx3 = x6;
dx4 = 2*n*x5+3*x1*n^2;
dx5 = -2*n*x4;
dx6 = -x3*n^2;

Pp = F*p+p*F'+Q;

dx = [dx1; dx2; dx3; dx4; dx5; dx6];
dp = [Pp(:,1); Pp(:,2); Pp(:,3); Pp(:,4); Pp(:,5); Pp(:,6)];
Xp = [dx; dp];