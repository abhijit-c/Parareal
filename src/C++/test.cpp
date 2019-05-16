#include <iostream>
#include <Eigen/Dense>

#include <math.h>
#include <stdio.h>
#include <omp.h>

#include "ode_integrators/integrators.h"

typedef Eigen::VectorXd Evec;
typedef Eigen::MatrixXd Emat;

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

  time_stepper course;
  course.dt = .1;
  course.F = std::function<int(ode_system, double, Evec &)>(&forward_euler);
  course.F_steps = std::function<int(ode_system, double, Emat &)>(&forward_euler_steps);

  time_stepper fine;
  fine.dt = .1*.1;
  fine.F = std::function<int(ode_system, double, Evec &)>(&forward_euler);
  fine.F_steps = std::function<int(ode_system, double, Emat &)>(&forward_euler_steps);

  Emat yf(ode.num_steps(.1), ode.dimension); 
  double tt = omp_get_wtime();

  parareal(ode, course, fine, yf);

  tt = omp_get_wtime() - tt;
  int k = yf.rows();
  printf("err=%e, computed in %f\n", exp(4) - yf(k-1,0), tt);
}
