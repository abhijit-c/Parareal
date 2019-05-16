#include <omp.h>
#include <iostream>
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
  Eigen::MatrixXd yfine(csteps+1, D);
  Eigen::MatrixXd ycourse(csteps+1, D);
  yfine.row(0) = sys.y0; ycourse.row(0) = sys.y0; yf.row(0) = sys.y0;
  for (int k = 0; k < 1; k++)
  {
    for (int n = 0; n < csteps; n++)
    { //Compute yf(n) = fine(yf(n-1)) w/
      ode_system para = sys;
      para.t_init = sys.t_init + course.dt*n;
      para.t_final = sys.t_init + course.dt*(n+1);
      para.y0 = yf.row(n);
      fine.integrate(para, temp);
      yfine.row(n+1) = temp;
    }
    std::cout << "FINE\n" << yfine << "\nFINE" << std::endl;
    for (int n = 0; n < csteps; n++)
    { // Predict w/ course operator, correct with fine.
      ode_system para = sys;
      para.t_init = sys.t_init + course.dt*n;
      para.t_final = sys.t_init + course.dt*(n+1);
      para.y0 = yf.row(n);
      course.integrate(para, temp);
      ycourse.row(n+1) = temp;
      yf.row(n+1) = ycourse.row(n+1)+yfine.row(n+1)-prev_course.row(n+1);
      prev_course.row(n+1) = ycourse.row(n+1);
    }
    std::cout << "COURSE\n" << prev_course << "\nCOURSE" << std::endl;
  }
  return 0;
}
