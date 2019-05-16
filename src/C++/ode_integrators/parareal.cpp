#include <omp.h>
#include "integrators.h"

int parareal(ode_system sys, 
             time_stepper course, time_stepper fine, 
             Eigen::VectorXd &yf)
{
  /*
  int D = sys.dimension, csteps = sys.num_steps(course.dt),
                         fsteps = sys.num_steps(fine.dt);
  Eigen::MatrixXd y_c(sys.num_steps(course.dt), D); 
  course.s_integrate(sys, y_c);
  */
  return 0;
}
