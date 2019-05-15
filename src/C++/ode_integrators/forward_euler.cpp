#include "integrators.h"

int forward_euler(ode_system sys, double dt, Eigen::VectorXd &yf)
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
