function F = calcularMatrizF(mi, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
clc
clear all

syms muu x y z xdot ydot zdot r0 r0dot theta0 theta0dot

f8 = (r0*theta0dot^2)-muu/(r0^2);
f10 = -2*r0dot*theta0dot/r0;

f4 = 2*theta0dot*ydot+f10*y+x*(theta0dot^2)-(muu*(r0+x)/(((r0+x)^2+y^2+z^2)^(3/2))) +muu/(r0^2);
f5 = -2*theta0dot*xdot-f10*x+y*(theta0dot^2)-(muu*y/(((r0+x)^2+y^2+z^2)^(3/2)));
f6 = -(muu*z/(((r0+x)^2+y^2+z^2)^(3/2)));

df4_dx1 = diff(f4, x);
df4_dx2 = diff(f4, y);
df4_dx3 = diff(f4, z);
df4_dx4 = diff(f4, xdot);
df4_dx5 = diff(f4, ydot);
df4_dx6 = diff(f4, zdot);
df4_dx7 = diff(f4, r0);
df4_dx8 = diff(f4, r0dot);
df4_dx9 = diff(f4, theta0);
df4_dx10 = diff(f4, theta0dot);

df5_dx1 = diff(f5, x);
df5_dx2 = diff(f5, y);
df5_dx3 = diff(f5, z);
df5_dx4 = diff(f5, xdot);
df5_dx5 = diff(f5, ydot);
df5_dx6 = diff(f5, zdot);
df5_dx7 = diff(f5, r0);
df5_dx8 = diff(f5, r0dot);
df5_dx9 = diff(f5, theta0);
df5_dx10 = diff(f5, theta0dot);

df6_dx1 = diff(f6, x);
df6_dx2 = diff(f6, y);
df6_dx3 = diff(f6, z);
df6_dx4 = diff(f6, xdot);
df6_dx5 = diff(f6, ydot);
df6_dx6 = diff(f6, zdot);
df6_dx7 = diff(f6, r0);
df6_dx8 = diff(f6, r0dot);
df6_dx9 = diff(f6, theta0);
df6_dx10 = diff(f6, theta0dot);

df8_dx1 = diff(f8, x);
df8_dx2 = diff(f8, y);
df8_dx3 = diff(f8, z);
df8_dx4 = diff(f8, xdot);
df8_dx5 = diff(f8, ydot);
df8_dx6 = diff(f8, zdot);
df8_dx7 = diff(f8, r0);
df8_dx8 = diff(f8, r0dot);
df8_dx9 = diff(f8, theta0);
df8_dx10 = diff(f8, theta0dot);

df10_dx1 = diff(f10, x)
df10_dx2 = diff(f10, y)
df10_dx3 = diff(f10, z)
df10_dx4 = diff(f10, xdot)
df10_dx5 = diff(f10, ydot)
df10_dx6 = diff(f10, zdot)
df10_dx7 = diff(f10, r0)
df10_dx8 = diff(f10, r0dot)
df10_dx9 = diff(f10, theta0)
df10_dx10 = diff(f10, theta0dot)

F = [0 0 0 1 0 0 0 0 0 0;
     0 0 0 0 1 0 0 0 0 0;
     0 0 0 0 0 1 0 0 0 0;
     df4_dx1 df4_dx2 df4_dx3 df4_dx4 df4_dx5 df4_dx6 df4_dx7 df4_dx8 df4_dx9 df4_dx10;
     df5_dx1 df5_dx2 df5_dx3 df5_dx4 df5_dx5 df5_dx6 df5_dx7 df5_dx8 df5_dx9 df5_dx10;
     df6_dx1 df6_dx2 df6_dx3 df6_dx4 df6_dx5 df6_dx6 df6_dx7 df6_dx8 df6_dx9 df6_dx10;
     0 0 0 0 0 0 0 1 0 0;
     df8_dx1 df8_dx2 df8_dx3 df8_dx4 df8_dx5 df8_dx6 df8_dx7 df8_dx8 df8_dx9 df8_dx10;
     0 0 0 0 0 0 0 0 0 1;
     df10_dx1 df10_dx2 df10_dx3 df10_dx4 df10_dx5 df10_dx6 df10_dx7 df10_dx8 df10_dx9 df10_dx10];

% disp(F);
end



 