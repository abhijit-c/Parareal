#include "integrators.h"

int rk4(ode_system sys, double dt, double *yf)
{
  int D = sys.dimension;
  for (int k = 0; k < D; k++)
  {
    yf[k] = sys.y0[k];
  }

  double yt[D], dydt[4*D];
  double t = sys.t_init;
  while (t + dt/2 < sys.t_final)
  {
    //k_1 = f(t, y)
    sys.f(t, yf, dydt);
    //k_2 = f(t + dt/2, y + k1/2)
    for (int k = 0; k < D; k++) { yt[k] = yf[k] + dt*dydt[k]/2; }
    sys.f(t + dt/2, yt, dydt+D);
    //k_3 = f(t + dt/2, y + k2/2)
    for (int k = 0; k < D; k++) { yt[k] = yf[k] + dt*dydt[k+D]/2; }
    sys.f(t + dt/2, yt, dydt+2*D);
    //k_4 = f(t + dt, y + k3)
    for (int k = 0; k < D; k++) { yt[k] = yf[k] + dt*dydt[k+2*D]; }
    sys.f(t + dt, yt, dydt+3*D);

    for (int k = 0; k < D; k++)
    { //y_{n+1} = y_n + k_1/6 + k_2/3 + k_3/3 + k_4/6
      yf[k] = yf[k]+dt*((dydt[k]+dydt[k+3*D])/6 + (dydt[k+D] + dydt[k+2*D])/3);
    }
    t = t + dt;
  }

  return 0;
}
