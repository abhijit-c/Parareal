#include "integrators.h"

// RK2 a.k.a. the midpoint method.
int rk2(ode_system &sys, double dt, Eigen::VectorXd &yf)
{
  int D = sys.dimension;
  yf = sys.y0;

  double t = sys.t_init;
  Eigen::VectorXd yt(D), k1(D), k2(D);
  while (t + dt/2 < sys.t_final)
  {
    sys.f(t, yf, k1);
    yt = yf + 0.5*dt*k1; sys.f(t + 0.5*dt, yt, k2);

    yf = yf + dt*k2;
    t = t + dt;
  }
  return 0;
}

int rk2_allt(ode_system &sys, double dt, Eigen::MatrixXd &steps)
{
  int D = sys.dimension, step = 1;
  Eigen::VectorXd yf(D), yt(D), k1(D), k2(D);
  yf = sys.y0; steps.row(0) = sys.y0;

  double t = sys.t_init;
  while (t + dt/2 < sys.t_final)
  {
    sys.f(t, yf, k1);
    yt = yf + 0.5*dt*k1; sys.f(t + 0.5*dt, yt, k2);

    yf = yf + dt*k2; steps.row(step) = yf;
    t = t + dt; step = step + 1;
  }
  return 0;
}
