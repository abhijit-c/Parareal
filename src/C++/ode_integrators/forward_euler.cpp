#include <math.h>
#include "integrators.h"
#include <stdio.h>

int forward_euler(ode_system &sys, double dt, Eigen::VectorXd &yf)
{
  yf = sys.y0;
  Eigen::VectorXd dydt(sys.dimension);
  double t = sys.t_init;
  while (t + dt/2 < sys.t_final)
  {
    sys.f(t, yf, dydt);
    yf = yf + dt*dydt;
    t = t + dt;
  }
  return 0;
}

int forward_euler_allt(ode_system &sys, double dt, Eigen::MatrixXd &steps)
{
  Eigen::VectorXd yf(sys.dimension), dydt(sys.dimension);

  yf = sys.y0; steps.row(0) = sys.y0;
  // Temporal steps initialization
  int step = 1; double t = sys.t_init;
  while (t + dt/2 < sys.t_final)
  {
    sys.f(t, yf, dydt); yf = yf + dt*dydt;
    steps.row(step) = yf;
    t = t + dt; step = step + 1;
  }
  return 0;
}
