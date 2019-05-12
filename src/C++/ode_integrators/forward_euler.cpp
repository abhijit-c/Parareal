#include "integrators.h"

int forward_euler(ode_system sys, double dt, double *yf)
{
  for (int k = 0; k < sys.dimension; k++)
  {
    yf[k] = sys.y0[k];
  }

  double dydt[sys.dimension];
  double t = sys.t_init;
  while (t + dt/2 < sys.t_final)
  {
    sys.f(t, yf, dydt);
    for (int k = 0; k < sys.dimension; k++)
    {
      yf[k] = yf[k] + dt*dydt[k];
    }
    t = t + dt;
  }
  return 0;
}
