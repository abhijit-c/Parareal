#include <iostream>
#include <fstream>

#include <math.h>
#include <stdio.h>
#include <omp.h>

#include "ode_integrators/integrators.h"
#include "Eigen/Dense"

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
  int csteps = ode.num_steps(course.dt), fsteps = ode.num_steps(fine.dt);
  double true_sol = ode.y0(0)*exp(LAMBDA*ode.t_final), tt = 0.0;
  // Fine Solver
  Emat yf_fine(fsteps, ode.dimension);
  tt = omp_get_wtime();
  fine.integrate_allt(ode, yf_fine);
  tt = omp_get_wtime() - tt;
  //std::cout << yf_fine << std::endl;
  printf("Fine: %e, w/ err %e, in %fs\n", 
      yf_fine(fsteps-1,0), abs(true_sol-yf_fine(fsteps-1,0)), tt);

  // Parareal Solver
  Emat yf(csteps, ode.dimension);
  tt = omp_get_wtime();
  parareal(ode, course, fine, 4, yf);
  tt = omp_get_wtime() - tt;
  std::cout << yf << std::endl;
  printf("Parareal: %e, w/ err %e, in %fs\n", 
      yf(csteps-1,0), abs(true_sol-yf(csteps-1,0)), tt);

  //Pipelined Parareal Solver
  Emat yfp(csteps, ode.dimension);
  tt = omp_get_wtime();
  pipelined_parareal(ode, course, fine, 4, yfp);
  tt = omp_get_wtime() - tt;
  std::cout << yfp << std::endl;
  printf("Pipelined Parareal: %e, w/ err %e, in %fs\n", 
      yfp(csteps-1,0), abs(true_sol-yfp(csteps-1,0)), tt);
}

void write_parareal(ode_system &ode, time_stepper course, time_stepper fine)
{
  int csteps = ode.num_steps(course.dt);
  for (int k = 0; k < csteps - 1; k++)
  {
    printf("Starting Parareal-%d solve\n", k);
    Emat yf(csteps, ode.dimension);
    parareal(ode, course, fine, k, yf);
    FILE* fd = NULL;
    char filename[256];
    snprintf(filename, 256, "Graphical/Data/linear1d_fw%02d.out", k);
    fd = fopen(filename,"w+");
    if(NULL == fd) {
      printf("Error opening file \n");
      return;
    }
    for(int i = 0; i < csteps; ++i) { fprintf(fd, "%f\n", yf(i)); }
    fclose(fd);
  }

}

void write_fine(ode_system &ode, time_stepper course, time_stepper fine)
{
  int fsteps = ode.num_steps(fine.dt);
  printf("Beginning fine solve, %d steps\n", fsteps);
  Emat yf(fsteps, ode.dimension);
  fine.integrate_allt(ode, yf);
  FILE* fd = NULL;
  char filename[256];
  snprintf(filename, 256, "Graphical/Data/linear1d_fine.out");
  fd = fopen(filename,"w+");
  if(NULL == fd) {
    printf("Error opening file \n");
    return;
  }
  for(int i = 0; i < fsteps; ++i) { fprintf(fd, "%f\n", yf(i)); }
  fclose(fd);
}

int main(int argc, char **argv)
{
  ode_system ode;
  ode.dimension = 1; ode.t_init = 0; ode.t_final = 4;
  ode.y0 = Evec(1); ode.y0(0) = 1;
  ode.f = std::function<int(double, Evec&, Evec&)>(&rhs);


  time_stepper course; time_stepper fine;
  // Forward Euler 
  course.dt = .5;
  course.F = std::function<int(ode_system&, double, Evec &)>(&forward_euler);
  course.F_allt = std::function<int(ode_system&, double, Emat &)>(&forward_euler_allt);
  fine.dt = .00001;
  fine.F = std::function<int(ode_system&, double, Evec &)>(&forward_euler);
  fine.F_allt = std::function<int(ode_system&, double, Emat &)>(&forward_euler_allt);
  printf("Forward Euler tests:\n");
  //write_parareal(ode, course, fine);
  //write_fine(ode, course, fine);
  test_method(ode, course, fine);

  
  /*
  // rk4 Testing 
  course.dt = .5;
  course.F = std::function<int(ode_system&, double, Evec &)>(&rk4);
  course.F_allt = std::function<int(ode_system&, double, Emat &)>(&rk4_allt);
  fine.dt = .0000001;
  fine.F = std::function<int(ode_system&, double, Evec &)>(&rk4);
  fine.F_allt = std::function<int(ode_system&, double, Emat &)>(&rk4_allt);
  printf("RK4 tests:\n");
  test_method(ode, course, fine);
  */

  return 0;
}
