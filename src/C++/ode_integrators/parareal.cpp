#include <omp.h>
#include "integrators.h"

int parareal(ode_system sys, 
             time_stepper course, time_stepper fine, 
             Eigen::VectorXd &yf)
{
  int D = sys.dimension;
  Eigen::VectorXd y_c(D); 
  course.sstep_integrate(sys, y_c);
  return 0;
}
