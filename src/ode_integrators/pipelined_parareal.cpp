#include <omp.h>
#include "integrators.h"

int pipelined_parareal(ode_system &sys, time_stepper course, time_stepper fine, 
             int para_its, Eigen::MatrixXd &yf)
{
  int D = sys.dimension, csteps = sys.num_steps(course.dt);

  // Initialize omp locks, one for each processor.
  omp_lock_t lock[csteps-1];
  for (int i = 0; i < csteps - 1; i++) { omp_init_lock(&(lock[i])); }

  // Initialize course/fine temporary structures
  Eigen::MatrixXd ycourse(csteps, D), yfine(csteps, D), delta_y(csteps, D);
  ycourse.row(0) = sys.y0; yfine.row(0) = sys.y0; yf.row(0) = sys.y0;
  Eigen::VectorXd tt(csteps);
  for (int k = 0; k < csteps; k++) { tt(k) = sys.t_init + k*course.dt; }

  #pragma omp parallel
  { //Begin Pipelined Parareal.
    int N = omp_get_max_threads();

    #pragma omp for nowait
    for (int p = 0; p < N; p++)
    { // BEGIN initial course solve 
      ode_system temp_sys = sys;
      Eigen::VectorXd y_temp(D);
      if (p != 0)
      {
        temp_sys.t_init = 0; temp_sys.t_final = tt(p);
        course.integrate(temp_sys, y_temp);
        yf.row(p) = y_temp;
      }
      temp_sys.t_init = tt(p); temp_sys.t_final = tt(p+1);
      temp_sys.y0 = y_temp;
      course.integrate(temp_sys, y_temp);
      ycourse.row(p+1) = y_temp;
    } // END initial initial solve NOWAIT

    for (int k = 1; k < para_its; k++)
    { // BEGIN parareal iterations
      #pragma omp for nowait ordered // ORDERED IMPORTANT
      for (int p = 0; p < N; p++)
      { //BEGIN processor p computation
        ode_system temp_sys = sys;
        Eigen::VectorXd y_temp(D);

        // Compute fine solution and corrector term for pth node.
        temp_sys.t_init = tt(p); temp_sys.t_final = tt(p+1);
        temp_sys.y0 = yf.row(p);
        omp_set_lock(&(lock[p]));
        fine.integrate(temp_sys, y_temp);
        delta_y.row(p+1) = y_temp - ycourse.row(p+1);
        omp_unset_lock(&(lock[p]));

        #pragma omp ordered
        { // BEGIN ordered region
          temp_sys.t_init = tt(p); temp_sys.t_final = tt(p+1);
          temp_sys.y0 = yf.row(p);
          course.integrate(temp_sys, y_temp);
          ycourse.row(p+1) = y_temp;
          if (p != N-1)
          { // Corrector step.
            omp_set_lock(&(lock[p+1]));
            yf.row(p+1) = ycourse.row(p+1) + delta_y.row(p+1);
            omp_unset_lock(&(lock[p+1]));
          }
        } // END ordered region
      } // END processor p computation
    } // END parareal iterations
  } // END pipelined parareal 

  // Clean up space.
  for (int i = 0; i < csteps - 1; i++) { omp_destroy_lock(&(lock[i])); }

  return 0;
}
