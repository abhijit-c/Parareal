#include <omp.h>
#include "integrators.h"

int parareal(ode_system sys, 
             time_stepper course, time_stepper fine, 
             Eigen::VectorXd &yf)
{
  //int P = omp_get_max_threads();
  return 0;
}
