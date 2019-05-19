#include "integrators.h"
#include <omp.h>

int parareal(ode_system &sys, time_stepper course, time_stepper fine, 
             int para_its, Eigen::MatrixXd &yf)
{
  int D = sys.dimension, csteps = sys.num_steps(course.dt);
  /* Serially compute the course solution to sys */
  course.integrate_allt(sys, yf);
  /* Initialize containers for parareal */
  Eigen::MatrixXd ycourse = yf;
  Eigen::MatrixXd yfine(csteps, D); yfine.row(0) = sys.y0;
  for (int k = 0; k < para_its; k++)
  { // Begin Parareal Steps
    /* In parallel compute the fine iterates on top of the serial steps */
    #pragma omp parallel for
    for (int n = 0; n < csteps-1; n++)
    { 
      /* Construct fine ODE */
      ode_system para = sys;
      para.t_init = sys.t_init + course.dt*n;
      para.t_final = sys.t_init + course.dt*(n+1);
      para.y0 = yf.row(n);
      /* Solve and update yfine */
      Eigen::VectorXd temp;
      fine.integrate(para, temp);
      yfine.row(n+1) = temp;
    }
    for (int n = 0; n < csteps-1; n++)
    { // Predict w/ course operator, correct with fine.
      /* Construct predictor ODE system */
      ode_system para = sys;
      para.t_init = sys.t_init + course.dt*n;
      para.t_final = sys.t_init + course.dt*(n+1);
      para.y0 = yf.row(n);
      /* Correct with parareal iterative scheme */
      Eigen::VectorXd temp(D);
      course.integrate(para, temp);
      temp = temp.transpose() + yfine.row(n+1) - ycourse.row(n+1);
      yf.row(n+1) = temp;
      ycourse.row(n+1) = temp;
    }
  }
  return 0;
}
