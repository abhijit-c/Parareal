#include <iostream>
#include <fstream>

#include <math.h>
#include <stdio.h>
#include <omp.h>

#include "ode_integrators/integrators.h"
#include "Eigen/Dense"

typedef Eigen::VectorXd Evec;
typedef Eigen::MatrixXd Emat;

#define SIGMA 10
#define BETA 2.6667
#define RHO 28

// u' = (LAMBDA)*u
int rhs(double t, Evec &yf, Evec &dydt)
{
  dydt(0) = SIGMA*( yf(1)-yf(0) );
  dydt(1) = yf(0)*( RHO - yf(2) ) - yf(1);
  dydt(2) = yf(0)*yf(1) - BETA*yf(2);
  return 0;
}

void test_method(ode_system &ode, time_stepper course, time_stepper fine)
{
  int csteps = ode.num_steps(course.dt), fsteps = ode.num_steps(fine.dt);
  double tt;
  // Fine Solver
  Emat yf_fine(fsteps, ode.dimension);
  tt = omp_get_wtime();
  fine.integrate_allt(ode, yf_fine);
  tt = omp_get_wtime() - tt;
  //printf("Fine Solution Computed in %fs\n", tt);
  //std::cout << yf_fine << std::endl;

  // Parareal Solver
  Emat yf(csteps, ode.dimension);
  tt = omp_get_wtime();
  parareal(ode, course, fine, 4, yf);
  tt = omp_get_wtime() - tt;
  //printf("Parareal Solution Computed in %fs\n", tt);
  std::cout << yf << std::endl;
  
  // Parareal Solver
  //tt = omp_get_wtime();
  //pipelined_parareal(ode, course, fine, 4, yf);
  //tt = omp_get_wtime() - tt;
  //printf("Parareal Solution Computed in %fs\n", tt);
  //std::cout << yf << std::endl;
}

int main(int argc, char **argv)
{
  ode_system ode;
  ode.dimension = 3; ode.t_init = 0; ode.t_final = 100;
  Evec y0(3); y0 << 1, 1, 1; ode.y0 = y0; 
  ode.f = std::function<int(double, Evec&, Evec&)>(&rhs);

  time_stepper course; time_stepper fine;
  /* rk4 Testing */
  course.dt = .05;
  course.F = std::function<int(ode_system&, double, Evec &)>(&rk4);
  course.F_allt = std::function<int(ode_system&, double, Emat &)>(&rk4_allt);
  fine.dt = .0001;
  fine.F = std::function<int(ode_system&, double, Evec &)>(&rk4);
  fine.F_allt = std::function<int(ode_system&, double, Emat &)>(&rk4_allt);
  printf("RK4 tests:\n");
  test_method(ode, course, fine);

  return 0;
}
