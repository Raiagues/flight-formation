
syms mu x y z xdot ydot zdot r0 r0dot theta0 theta0dot

dx7 = r0dot;                  
dx8 = r0*theta0dot^2 - mu/(r0^2);  
dx9 = theta0dot;                  
dx10 = -2*r0dot*theta0dot/r0;       

dx1 = xdot;
dx2 = ydot;
dx3 = zdot;
dx4 = 2*theta0dot*ydot + dx10*y + (theta0dot^2)*x - mu*(r0+x)/((r0+x)^2+y^2+z^2)^(3/2) + mu/(r0^2);
dx5 = -2*theta0dot*xdot - dx10*x + (theta0dot^2)*y - mu*y/((r0+x)^2+y^2+z^2)^(3/2);
dx6 = -mu*z/((r0+x)^2+y^2+z^2)^(3/2);


df4_dx1 = diff(dx4, x);
df4_dx2 = diff(dx4, y);
df4_dx3 = diff(dx4, z);
df4_dx4 = diff(dx4, xdot);
df4_dx5 = diff(dx4, ydot);
df4_dx6 = diff(dx4, zdot);
df4_dx7 = diff(dx4, r0);
df4_dx8 = diff(dx4, r0dot);
df4_dx9 = diff(dx4, theta0);
df4_dx10 = diff(dx4, theta0dot);

df5_dx1 = diff(dx5, x);
df5_dx2 = diff(dx5, y);
df5_dx3 = diff(dx5, z);
df5_dx4 = diff(dx5, xdot);
df5_dx5 = diff(dx5, ydot);
df5_dx6 = diff(dx5, zdot);
df5_dx7 = diff(dx5, r0);
df5_dx8 = diff(dx5, r0dot);
df5_dx9 = diff(dx5, theta0);
df5_dx10 = diff(dx5, theta0dot);

df6_dx1 = diff(dx6, x);
df6_dx2 = diff(dx6, y);
df6_dx3 = diff(dx6, z);
df6_dx4 = diff(dx6, xdot);
df6_dx5 = diff(dx6, ydot);
df6_dx6 = diff(dx6, zdot);
df6_dx7 = diff(dx6, r0);
df6_dx8 = diff(dx6, r0dot);
df6_dx9 = diff(dx6, theta0);
df6_dx10 = diff(dx6, theta0dot);

df8_dx1 = diff(dx8, x);
df8_dx2 = diff(dx8, y);
df8_dx3 = diff(dx8, z);
df8_dx4 = diff(dx8, xdot);
df8_dx5 = diff(dx8, ydot);
df8_dx6 = diff(dx8, zdot);
df8_dx7 = diff(dx8, r0);
df8_dx8 = diff(dx8, r0dot);
df8_dx9 = diff(dx8, theta0);
df8_dx10 = diff(dx8, theta0dot);

df10_dx1 = diff(dx10, x);
df10_dx2 = diff(dx10, y);
df10_dx3 = diff(dx10, z);
df10_dx4 = diff(dx10, xdot);
df10_dx5 = diff(dx10, ydot);
df10_dx6 = diff(dx10, zdot);
df10_dx7 = diff(dx10, r0);
df10_dx8 = diff(dx10, r0dot);
df10_dx9 = diff(dx10, theta0);
df10_dx10 = diff(dx10, theta0dot);


test = [0 0 0 1 0 0 0 0 0 0;
        0 0 0 0 1 0 0 0 0 0;
        0 0 0 0 0 1 0 0 0 0;
        df4_dx1 df4_dx2 df4_dx3 df4_dx4 df4_dx5 df4_dx6 df4_dx7 df4_dx8 df4_dx9 df4_dx10;
        df5_dx1 df5_dx2 df5_dx3 df5_dx4 df5_dx5 df5_dx6 df5_dx7 df5_dx8 df5_dx9 df5_dx10;
        df6_dx1 df6_dx2 df6_dx3 df6_dx4 df6_dx5 df6_dx6 df6_dx7 df6_dx8 df6_dx9 df6_dx10;
        0 0 0 0 0 0 0 1 0 0;
        df8_dx1 df8_dx2 df8_dx3 df8_dx4 df8_dx5 df8_dx6 df8_dx7 df8_dx8 df8_dx9 df8_dx10;
        0 0 0 0 0 0 0 0 0 1;
        df10_dx1 df10_dx2 df10_dx3 df10_dx4 df10_dx5 df10_dx6 df10_dx7 df10_dx8 df10_dx9 df10_dx10]
 
disp(test)
