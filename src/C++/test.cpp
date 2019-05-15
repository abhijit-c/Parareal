#include <Eigen/Dense>
#include <math.h>
#include <stdio.h>
#include <omp.h>
#include "ode_integrators/integrators.h"

typedef Eigen::VectorXd Evec;

int rhs(double t, Evec &yf, Evec &dydt)
{
  dydt(0) = yf(0);
  return 0;
}

int main()
{
  ode_system ode;
  ode.dimension = 1; ode.t_init = 0; ode.t_final = 4;
  ode.y0 = Evec(1); ode.y0(0) = 1;
  ode.f = std::function<int(double, Evec&, Evec&)>(&rhs);

  Evec yf(1); yf(0) = 0;

  forward_euler(ode, .001, yf);
  printf("FW: Computed: %f, error = %e\n", yf[0], exp(ode.t_final)-yf[0]);

  rk4(ode, .001, yf);
  printf("RK4: Computed: %f, error = %e\n", yf[0], exp(ode.t_final)-yf[0]);
}
