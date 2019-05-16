#include <omp.h>
#include "integrators.h"

int parareal(ode_system sys, 
             time_stepper course, time_stepper fine, 
             Eigen::MatrixXd &yf)
{
  int P = omp_get_max_threads();
  int D = sys.dimension, csteps = sys.num_steps(course.dt),
                         fsteps = sys.num_steps(fine.dt);
  course.s_integrate(sys, yf);
  Eigen::VectorXd temp;
  Eigen::MatrixXd prev_course = yf;
  Eigen::MatrixXd yfine(csteps, D);
  Eigen::MatrixXd ycourse(csteps, D);
  yfine.row(0) = sys.y0; ycourse.row(0) = sys.y0; yf.row(0) = sys.y0;
  for (int k = 0; k < P; k++)
  {
    for (int n = 0; n < csteps-1; n++)
    { //Compute yf(n) = fine(yf(n-1)) w/
      ode_system para = sys;
      para.t_init = para.t_init + course.dt*n;
      para.t_final = para.t_final + course.dt*(n+1);
      para.y0 = yf.row(n);
      fine.integrate(para, temp);
      yfine.row(n+1) = temp;
    }
    for (int n = 0; n < csteps-1; n++)
    { //Compute yf(n) = fine(yf(n-1)) w/
      ode_system para = sys;
      para.t_init = para.t_init + course.dt*n;
      para.t_final = para.t_final + course.dt*(n+1);
      para.y0 = ycourse.row(n);
      course.integrate(para, temp);
      ycourse.row(n+1) = temp;
      yf.row(n+1) = ycourse.row(n+1)+yfine.row(n+1)-prev_course.row(n+1);
    }
    prev_course = ycourse;
  }
  return 0;
}
