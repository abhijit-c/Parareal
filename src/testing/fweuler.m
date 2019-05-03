function [t_vals, y_vals] = fweuler(ode, y0, t_range, dt)
% fweuler is a function which computes an approximate solution to the first
% order ODE system y' = ode(t,y) with initial condition y0 using the foward
% euler scheme:
% \[ y(t_{n+1}) = y(t_n) + dt*ode(y_n + y(t_n)) \]
% The solver is agnositic to both the dimension of system and data points.
% ------------------------------------------------------------------------------
% Parameter List:
% INPUT:
% ode    : Right hand side of ODE, takes input (t,y). If [m,n] = size(y), m
%          describes the number of equations in the system, and n descirbes the
%          dimensionality of the points, i.e. n = 2 => y_i is in \mathbb{R}^2
% y0     : Initial condition at time t_range(0). This should have m rows
%          corresponding to the components of the system,  with columns
%          corresponding to the dimension
% t_range: Vector of size 2 describing time domain [initial_time, final_time]
% dt     : Desired time step
% OUTPUT:
% t_vals: The array [t_range(1):dt:t_range(2)].
% y_vals: The approximated solutions of the ode at each of the time values.
% ------------------------------------------------------------------------------
% Tested on GNU Octave 5.1.0. Matlab untested (should be fine).
% Author: Abhijit Chowdhary, Undergraduate @ New York University.
% Date  : 2019/04/19 (YYYY/MM/DD ISO 8601)
  [M,N] = size(y0);
  t_init = t_range(1); t_last = t_range(2);
  t_vals = [t_init:dt:t_last]';
  y = y0; y_vals = zeros(M*length(t_vals),N); y_vals(1:M, :) = y0;

  % Forward Euler temporal steps
  for k=1:length(t_vals)-1
    y = y + dt*( ode(t_vals(k),y) );
    y_vals(k*M+1:(k+1)*M, :) = y;
  end
end
