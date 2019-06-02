#include <stdio.h>
#include <omp.h>

#include <Eigen/Dense>
#include <Parareal/core.h>
#include <Parareal/parareal.h>
#include <Parareal/forward_euler.h>

#include "rhs_linear1d.h"

typedef Eigen::VectorXd Evec;
typedef Eigen::MatrixXd Emat;

int main(int argc, char **argv)
{
  int P = omp_get_max_threads();
  ode_system ode;
  ode.dimension = 1; ode.t_init = 0; ode.t_final = 0.5*P;
  ode.y0 = Evec(1); ode.y0(0) = 1;
  ode.f = std::function<int(double, Evec&, Evec&)>(&linear1d);

  printf("Using %d processors, testing strong scaling. \n", P);
  printf("Threads\tTime\n");
  for (int k = 2; k < P; k++)
  {
    omp_set_num_threads(k);

    time_stepper course; time_stepper fine;

    // Forward Euler 
    course.dt = 0.5;
    course.F = std::function<int(ode_system&, double, Evec &)>(&forward_euler);
    course.F_allt = std::function<int(ode_system&, double, Emat &)>(&forward_euler_allt);
    fine.dt = .0000001;
    fine.F = std::function<int(ode_system&, double, Evec &)>(&forward_euler);
    fine.F_allt = std::function<int(ode_system&, double, Emat &)>(&forward_euler_allt);

    // Test Parareal's Strong Scaling

    int csteps = ode.num_steps(course.dt);
    double tt = 0.0;

    // Parareal Solver
    Emat yf(csteps, ode.dimension);
    tt = omp_get_wtime();
    parareal(ode, course, fine, 6, yf);
    tt = omp_get_wtime() - tt;

    printf("%d\t%f\n", k, tt);
  }
  return 0;
}
