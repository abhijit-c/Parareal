#include <math.h>
#include <stdio.h>
#include <omp.h>

#include <iostream>

#include <Eigen/Dense>
#include <Parareal/core.h>
#include <Parareal/parareal.h>
#include <Parareal/forward_euler.h>

typedef Eigen::VectorXd Evec;
typedef Eigen::MatrixXd Emat;

#define LAMBDA 1

// u' = (LAMBDA)*u
int rhs(double t, Evec &yf, Evec &dydt)
{
  dydt(0) = LAMBDA*yf(0);
  return 0;
}

void test_method(ode_system &ode, time_stepper course, time_stepper fine)
{
  int csteps = ode.num_steps(course.dt);
  double tt = 0.0;
  // Parareal Solver
  Emat yf(csteps, ode.dimension);
  tt = omp_get_wtime();
  parareal(ode, course, fine, 4, yf);
  tt = omp_get_wtime() - tt;
  std::cout << yf << std::endl;
  printf("Time taken %f\n", tt);
}

int main(int argc, char **argv)
{
  int P = omp_get_max_threads();
  ode_system ode;
  ode.dimension = 1; ode.t_init = 0; ode.t_final = 0.5*P;
  ode.y0 = Evec(1); ode.y0(0) = 1;
  ode.f = std::function<int(double, Evec&, Evec&)>(&rhs);

  time_stepper course; time_stepper fine;

  // Forward Euler 
  course.dt = 0.5;
  course.F = std::function<int(ode_system&, double, Evec &)>(&forward_euler);
  course.F_allt = std::function<int(ode_system&, double, Emat &)>(&forward_euler_allt);
  fine.dt = .0000001;
  fine.F = std::function<int(ode_system&, double, Evec &)>(&forward_euler);
  fine.F_allt = std::function<int(ode_system&, double, Emat &)>(&forward_euler_allt);

  printf("Using %d processors to solve up to T = %f\n", P, ode.t_final);
  test_method(ode, course, fine);

  return 0;
}
