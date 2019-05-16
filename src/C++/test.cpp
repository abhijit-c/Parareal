#include <iostream>
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

  double tt = omp_get_wtime();
  forward_euler(ode, .001, yf);
  tt = omp_get_wtime() - tt;
  printf("FW: Computed: %f, error = %e\n, in %f s\n", 
      yf(0), exp(ode.t_final)-yf(0), tt);

  Eigen::MatrixXd M ( ode.num_steps(.001), ode.dimension );
  tt = omp_get_wtime(); 
  forward_euler_steps(ode, .001, M);
  tt = omp_get_wtime() - tt;
  int k = M.rows();
  printf("FW: Computed: %f, error = %e, in %f s\n", yf(0), 
      exp(ode.t_final)-yf(0), tt);

  tt = omp_get_wtime(); 
  rk4(ode, .001, yf);
  tt = omp_get_wtime() - tt;
  printf("RK4: Computed: %f, error = %e\n, in $f s\n", 
      yf[0], exp(ode.t_final)-yf(0), tt);
}
