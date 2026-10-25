function F = calcularMatrizF(mi, x1, x2, x3, x4, x5, x6)
clear all

syms mi x y z xdot ydot zdot n
% 
% mi = 3.986004418e14;
% Earth_r = 6378140;
% a = Earth_r + 600e3;
% n = sqrt(mi/a^3);

f4 = 2*n*ydot+3*x*n^2;
f5 = -2*n*xdot;
f6 = -z*n^2;

df4_dx1 = diff(f4, x);
df4_dx2 = diff(f4, y);
df4_dx3 = diff(f4, z);
df4_dx4 = diff(f4, xdot);
df4_dx5 = diff(f4, ydot);
df4_dx6 = diff(f4, zdot);

df5_dx1 = diff(f5, x);
df5_dx2 = diff(f5, y);
df5_dx3 = diff(f5, z);
df5_dx4 = diff(f5, xdot);
df5_dx5 = diff(f5, ydot);
df5_dx6 = diff(f5, zdot);

df6_dx1 = diff(f6, x);
df6_dx2 = diff(f6, y);
df6_dx3 = diff(f6, z);
df6_dx4 = diff(f6, xdot);
df6_dx5 = diff(f6, ydot);
df6_dx6 = diff(f6, zdot);


matriz = [0 0 0 1 0 0;
          0 0 0 0 1 0;
          0 0 0 0 0 1;
          df4_dx1 df4_dx2 df4_dx3 df4_dx4 df4_dx5 df4_dx6;
          df5_dx1 df5_dx2 df5_dx3 df5_dx4 df5_dx5 df5_dx6;
          df6_dx1 df6_dx2 df6_dx3 df6_dx4 df6_dx5 df6_dx6];
 


F = double(matriz);


end

