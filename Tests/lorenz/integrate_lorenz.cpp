#include <stdio.h>
#include <omp.h>

#include <iostream>

#include <Eigen/Dense>
#include <Parareal/core.h>
#include <Parareal/backward_euler.h>

#include "rhs_lorenz.h"

typedef Eigen::VectorXd Evec;
typedef Eigen::MatrixXd Emat;

int main(int argc, char **argv)
{
  ode_system ode;
  ode.dimension = 3; ode.t_init = 0; ode.t_final = 10;
  ode.y0 = Evec(1); ode.y0(0) = 1; ode.y0(1) = 0; ode.y0(2) = 0;
  ode.f = std::function<int(double, Evec&, Evec&)>(&lorenz);
  ode.J = std::function<int(double, Evec&, Emat&)>(&lorenz_jac);

  time_stepper solver; 

  // Backward Euler 
  solver.dt = 0.01;
  solver.F = std::function<int(ode_system&, double, Evec &)>(&backward_euler);
  solver.F_allt = std::function<int(ode_system&, double, Emat &)>(&backward_euler_allt);

  int steps = ode.num_steps(solver.dt);
  double tt = 0.0;

  // Parareal Solver
  Emat yf(steps, ode.dimension);
  tt = omp_get_wtime();
  solver.integrate_allt(ode, yf);
  tt = omp_get_wtime() - tt;

  std::cout << yf << std::endl;
  printf("Time taken %f\n", tt);

  return 0;
}
