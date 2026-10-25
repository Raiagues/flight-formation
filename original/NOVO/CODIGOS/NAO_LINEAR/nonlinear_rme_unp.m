function nrmv = nonlinear_rme_unp(t,nrmv0)

% RELATIVE MOTION EQUATIONS OF THE UNPERTURBED CASE

% Nonlinear relative motion variables (nrmv):
%   nrmv(1) = x coordinate (m)
%   nrmv(2) = y coordinate (m)
%   nrmv(3) = z coordinate (m)
%   nrmv(4) = x_dot coordinate (m/s)
%   nrmv(5) = y_dot coordinate (m/s)
%   nrmv(6) = z_dot coordinate (m/s)
%   nrmv(7) = Leader's radius in the ECI frame (m)
%   nrmv(8) = Leader's velocity in the ECI frame (m/s)
%   nrmv(9) = Leader's true anomaly (rad)
%   nrmv(10) = Leader's latitude argument rate (rad/s)

% Initial condition
global mi

x = nrmv0(1);
y = nrmv0(2);
z = nrmv0(3);

x_dot = nrmv0(4);
y_dot = nrmv0(5);
z_dot = nrmv0(6);

r0 = nrmv0(7);
r0_dot = nrmv0(8);

teta0 = nrmv0(9);
teta0_dot = nrmv0(10);

% Leader's acceleration in the ECI frame (m/s^2)
r0_2dot = r0*teta0_dot^2 - (mi/r0^2);

% Leader's true anomaly variation (rad/s^2)
teta0_2dot = -(2*r0_dot*teta0_dot)/r0;

% Relative accelerations (m/s^2)
x_2dot = 2*teta0_dot*y_dot + teta0_2dot*y + (teta0_dot^2)*x ...
          - ((mi*(r0 + x))/((r0 + x)^2 + y^2 + z^2)^(3/2)) + (mi/r0^2);
y_2dot = -2*teta0_dot*x_dot - teta0_2dot*x + (teta0_dot^2)*y ...
          - ((mi*y)/((r0 + x)^2 + y^2 + z^2)^(3/2));  
z_2dot = - ((mi*z)/((r0 + x)^2 + y^2 + z^2)^(3/2));

% Output
nrmv = [x_dot; y_dot; z_dot; x_2dot; y_2dot; z_2dot; r0_dot; r0_2dot; teta0_dot; teta0_2dot];
 
end

