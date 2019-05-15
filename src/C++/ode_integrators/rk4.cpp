#include "integrators.h"

int rk4(ode_system sys, double dt, Eigen::VectorXd &yf)
{
  int D = sys.dimension;
  yf = sys.y0;

  double t = sys.t_init;
  Eigen::VectorXd yt(D), k1(D), k2(D), k3(D), k4(D);
  while (t + dt/2 < sys.t_final)
  {
    sys.f(t, yf, k1);
    yt = yf + 0.5*dt*k1; sys.f(t + 0.5*dt, yt, k2);
    yt = yf + 0.5*dt*k2; sys.f(t + 0.5*dt, yt, k3);
    yt = yf + dt*k3;     sys.f(t + dt,     yt, k4);

    yf = yf + dt*( k1/6 + k2/3 + k3/3 + k4/6 );
    t = t + dt;
  }

  return 0;
}
